#!/usr/bin/env python3
import argparse
import math
import os
import sys
import numpy as np

DPI = 100
CBAR_PX = 120   # pixels reserved for colorbar panel
SCROLLBAR_W = 18  # scrollbar widget thickness
QUIT_H = 34       # quit button row height
DECO_H = 190      # OS title bar + taskbar allowance


def getScreenSize():
    """Return (width_px, height_px) via tkinter; fall back to 1920×1080."""
    try:
        import tkinter as tk
        root = tk.Tk()
        root.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        root.destroy()
        return w, h
    except Exception:
        return 1920, 1080


def blockAverage(arr, factor):
    """Block-average arr in factor×factor tiles using nanmean."""
    ny, nx = arr.shape[:2]
    ny2 = (ny // factor) * factor
    nx2 = (nx // factor) * factor
    a = arr[:ny2, :nx2]
    if a.ndim == 2:
        return np.nanmean(
            a.reshape(ny2 // factor, factor, nx2 // factor, factor),
            axis=(1, 3))
    nb = a.shape[2]
    return np.nanmean(
        a.reshape(ny2 // factor, factor, nx2 // factor, factor, nb),
        axis=(1, 3))


def readBand(ds, b):
    band = ds.GetRasterBand(b)
    data = band.ReadAsArray().astype(np.float32)
    nd = band.GetNoDataValue()
    if nd is not None:
        data[data == nd] = np.nan
    return data


def extractProfile(dec, r0, c0, r1, c1):
    """Sample dec along the line from (r0,c0) to (r1,c1) using bilinear interpolation."""
    length = max(int(np.hypot(r1 - r0, c1 - c0)), 1) + 1
    rows = np.linspace(r0, r1, length)
    cols = np.linspace(c0, c1, length)
    try:
        from scipy.ndimage import map_coordinates
        if dec.ndim == 2:
            vals = map_coordinates(np.nan_to_num(dec), [rows, cols], order=1, cval=np.nan)
        else:
            vals = np.stack(
                [map_coordinates(np.nan_to_num(dec[:, :, i]), [rows, cols],
                                 order=1, cval=np.nan)
                 for i in range(dec.shape[2])], axis=-1)
    except ImportError:
        # nearest-neighbour fallback
        ri = np.clip(np.round(rows).astype(int), 0, dec.shape[0] - 1)
        ci = np.clip(np.round(cols).astype(int), 0, dec.shape[1] - 1)
        vals = dec[ri, ci] if dec.ndim == 2 else dec[ri, ci, :]
    dist = np.linspace(0, np.hypot(c1 - c0, r1 - r0), length)
    return dist, vals


def openProfileWindow(dist, vals_or_list, p0, p1, titles=None, parent=None, pos=None):
    """Show profile(s) in a Toplevel window.

    vals_or_list: single ndarray or list of ndarrays (one per image).
    titles: optional list of per-subplot titles.
    """
    import tkinter as tk
    from tkinter import ttk
    import matplotlib.figure as mfig
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    vals_list = vals_or_list if isinstance(vals_or_list, list) else [vals_or_list]
    n = len(vals_list)
    if titles is None:
        titles = [None] * n

    WIN_W = 800
    WIN_H = 200 + 200 * n

    win = tk.Toplevel()
    win.title(f'Profile  ({p0[1]}, {p0[0]}) → ({p1[1]}, {p1[0]})')

    if pos is not None:
        win.geometry(f'{WIN_W}x{WIN_H}+{pos[0]}+{pos[1]}')
    elif parent is not None:
        parent.update_idletasks()
        px = parent.winfo_x()
        pw = parent.winfo_width()
        sw2 = parent.winfo_screenwidth()
        sh2 = parent.winfo_screenheight()
        x = max(0, min(px + pw + 10, sw2 - WIN_W))
        y = max(0, (sh2 - WIN_H) // 2)
        win.geometry(f'{WIN_W}x{WIN_H}+{x}+{y}')
    else:
        win.geometry(f'{WIN_W}x{WIN_H}')

    fig = mfig.Figure(figsize=(8, WIN_H / DPI), dpi=DPI)
    for i, (vals, ttl) in enumerate(zip(vals_list, titles)):
        ax = fig.add_subplot(n, 1, i + 1)
        if vals.ndim == 1:
            ax.plot(dist, vals, color='steelblue')
        else:
            for j in range(vals.shape[1]):
                ax.plot(dist, vals[:, j], label=f'Band {j + 1}')
            ax.legend(fontsize=7)
        ax.set_xlabel('Distance (pixels)')
        ax.set_ylabel('Value')
        hdr = f'{ttl}: ' if ttl else ''
        ax.set_title(f'{hdr}col={p0[1]}, row={p0[0]}  →  col={p1[1]}, row={p1[0]}'
                     f'   ({dist[-1]:.1f} px)', fontsize=8)
        ax.grid(True, alpha=0.4)
    fig.tight_layout()

    ttk.Button(win, text='Close', command=win.destroy).pack(side='bottom', pady=4)
    mpl = FigureCanvasTkAgg(fig, master=win)
    mpl.draw()
    mpl.get_tk_widget().pack(fill='both', expand=True)


# -----------------------------------------------------------------------
# Non-scroll path: single matplotlib figure with embedded colorbar
# -----------------------------------------------------------------------

def makeFigure(dec, title, cmap, vmin, vmax, is_rgb):
    """Build a matplotlib Figure sized exactly to the image in pixels."""
    import matplotlib.figure as mfig

    ny, nx = dec.shape[:2]
    fig_w_px = nx if is_rgb else nx + CBAR_PX
    fig = mfig.Figure(figsize=(fig_w_px / DPI, ny / DPI), dpi=DPI)

    if is_rgb:
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(dec, interpolation='nearest', aspect='equal')
    else:
        cbar_frac = CBAR_PX / fig_w_px
        ax_right = 1 - cbar_frac - 0.02
        ax = fig.add_axes([0, 0, ax_right, 1])
        cax = fig.add_axes([ax_right + 0.03, 0.05, 0.04, 0.9])
        im = ax.imshow(dec, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation='nearest', aspect='equal')
        fig.colorbar(im, cax=cax)

    ax.set_title(title, fontsize=8)
    ax.axis('off')
    return fig, fig_w_px, ny


# -----------------------------------------------------------------------
# Scroll path: PhotoImage for fast panning + separate colorbar figure
# -----------------------------------------------------------------------

def decToPhoto(dec, cmap, vmin, vmax, is_rgb):
    """Render dec to a PIL PhotoImage using the matplotlib colormap."""
    from PIL import Image, ImageTk
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors

    if is_rgb:
        arr = (np.clip(dec, 0, 1) * 255).astype(np.uint8)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        rgba = mcm.get_cmap(cmap)(norm(np.nan_to_num(dec, nan=vmin)))
        arr = (rgba[:, :, :3] * 255).astype(np.uint8)

    return ImageTk.PhotoImage(Image.fromarray(arr))


def makeColorbarFig(cmap, vmin, vmax, height_px):
    """Standalone colorbar figure for the scroll layout."""
    import matplotlib.figure as mfig
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors

    fig = mfig.Figure(figsize=(CBAR_PX / DPI, height_px / DPI), dpi=DPI)
    cax = fig.add_axes([0.25, 0.05, 0.35, 0.9])
    sm = mcm.ScalarMappable(cmap=mcm.get_cmap(cmap),
                             norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    fig.colorbar(sm, cax=cax)
    return fig


def bindScroll(tk_canvas):
    """Bind mouse-wheel scroll for Windows/Mac and Linux."""
    def _y(event):
        tk_canvas.yview_scroll(int(-1 * (event.delta / 120)), 'units')
    def _x(event):
        tk_canvas.xview_scroll(int(-1 * (event.delta / 120)), 'units')
    tk_canvas.bind('<MouseWheel>', _y)
    tk_canvas.bind('<Shift-MouseWheel>', _x)
    tk_canvas.bind('<Button-4>', lambda e: tk_canvas.yview_scroll(-1, 'units'))
    tk_canvas.bind('<Button-5>', lambda e: tk_canvas.yview_scroll(1, 'units'))
    tk_canvas.bind('<Shift-Button-4>', lambda e: tk_canvas.xview_scroll(-1, 'units'))
    tk_canvas.bind('<Shift-Button-5>', lambda e: tk_canvas.xview_scroll(1, 'units'))


# -----------------------------------------------------------------------
# Main display
# -----------------------------------------------------------------------

def showImage(image_defs, sw, sh):
    """Display 1–3 same-size images side by side with a floating control palette.

    image_defs: list of dicts, each with keys:
        dec (ndarray), title (str), cmap (str), vmin (float), vmax (float), is_rgb (bool)
    All images must share the same (ny, nx) shape.
    """
    import tkinter as tk
    from tkinter import ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    n_imgs = len(image_defs)
    ny, nx = image_defs[0]['dec'].shape[:2]

    root = tk.Tk()
    root.title(' | '.join(f'{i+1}) {os.path.basename(d["title"])}'
                          for i, d in enumerate(image_defs)))

    # ---- command palette (separate floating window) ----
    palette = tk.Toplevel(root)
    palette.title('Controls')
    palette.resizable(False, False)
    palette.protocol('WM_DELETE_WINDOW', root.destroy)

    pick_active         = [False]
    profile_active      = [False]
    col_active          = [False]
    row_active          = [False]
    lines_visible       = [True]
    overlay_set_visible = [None]
    profile_pts         = []
    lineplot_state      = {'win': None, 'axes': None, 'fig': None, 'canvas': None, 'mode': None}

    btn_col = ttk.Frame(palette)
    btn_col.pack(side='top', fill='x', padx=4, pady=(4, 2))

    pick_btn    = ttk.Button(btn_col, text='Pick')
    profile_btn = ttk.Button(btn_col, text='Profile')
    col_btn     = ttk.Button(btn_col, text='Col Plot')
    row_btn     = ttk.Button(btn_col, text='Row Plot')
    lines_btn   = ttk.Button(btn_col, text='Lines ✓')
    quit_btn    = ttk.Button(btn_col, text='Quit', command=root.destroy)
    for btn in (pick_btn, profile_btn, col_btn, row_btn, lines_btn, quit_btn):
        btn.pack(side='top', fill='x', pady=2, padx=2)

    status_var = tk.StringVar(value='Ready')
    status_lbl = ttk.Label(palette, textvariable=status_var, anchor='nw', wraplength=130)
    status_lbl.pack(side='bottom', fill='both', expand=True, padx=6, pady=(0, 4))

    def deactivate_all():
        pick_active[0] = False
        profile_active[0] = False
        col_active[0] = False
        row_active[0] = False
        pick_btn.config(text='Pick')
        profile_btn.config(text='Profile')
        col_btn.config(text='Col Plot')
        row_btn.config(text='Row Plot')

    def toggle_pick():
        if pick_active[0]:
            deactivate_all()
        else:
            deactivate_all()
            pick_active[0] = True
            pick_btn.config(text='Pick ●')
            profile_pts.clear()
    pick_btn.config(command=toggle_pick)

    def toggle_profile():
        if profile_active[0]:
            deactivate_all()
        else:
            deactivate_all()
            profile_active[0] = True
            profile_btn.config(text='Profile ●')
            profile_pts.clear()
            status_var.set('  Profile: click first point')
    profile_btn.config(command=toggle_profile)

    def toggle_col():
        if col_active[0]:
            deactivate_all()
        else:
            deactivate_all()
            col_active[0] = True
            col_btn.config(text='Col Plot ●')
            status_var.set('  Col Plot: click a pixel to plot that column')
    col_btn.config(command=toggle_col)

    def toggle_row():
        if row_active[0]:
            deactivate_all()
        else:
            deactivate_all()
            row_active[0] = True
            row_btn.config(text='Row Plot ●')
            status_var.set('  Row Plot: click a pixel to plot that row')
    row_btn.config(command=toggle_row)

    def toggle_lines():
        lines_visible[0] = not lines_visible[0]
        lines_btn.config(text='Lines ✓' if lines_visible[0] else 'Lines ✗')
        if overlay_set_visible[0]:
            overlay_set_visible[0](lines_visible[0])
    lines_btn.config(command=toggle_lines)

    def _lineplotWinAlive():
        win = lineplot_state['win']
        if win is None:
            return False
        try:
            return bool(win.winfo_exists())
        except Exception:
            return False

    def openOrReuseLineplot(mode):
        import matplotlib.figure as mfig
        if not _lineplotWinAlive() or lineplot_state['mode'] != mode:
            if _lineplotWinAlive():
                lineplot_state['win'].destroy()
            win = tk.Toplevel()
            win.title('Column Plots' if mode == 'col' else 'Row Plots')
            WIN_W = 800
            BTN_H = 46
            WIN_H = 200 + 200 * n_imgs
            x, y = nextPlotGeometry(WIN_W, WIN_H)
            win.geometry(f'{WIN_W}x{WIN_H}+{x}+{y}')
            fig = mfig.Figure(figsize=(8, (WIN_H - BTN_H) / DPI), dpi=DPI)
            axes = []
            for i, idef in enumerate(image_defs):
                ax = fig.add_subplot(n_imgs, 1, i + 1)
                ax.grid(True, alpha=0.4)
                ax.set_xlabel('Row index' if mode == 'col' else 'Column index')
                ax.set_ylabel('Value')
                ax.set_title(idef['title'], fontsize=8)
                axes.append(ax)
            fig.tight_layout()
            btn_frame = ttk.Frame(win)
            btn_frame.pack(side='bottom', fill='x', padx=4, pady=6)

            def save_plot():
                from tkinter import filedialog
                path = filedialog.asksaveasfilename(
                    parent=win,
                    defaultextension='.png',
                    filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'),
                               ('SVG', '*.svg'), ('All files', '*.*')])
                if path:
                    lineplot_state['fig'].savefig(path, bbox_inches='tight')

            ttk.Button(btn_frame, text='Save', command=save_plot).pack(side='left', padx=4)
            ttk.Button(btn_frame, text='Close', command=win.destroy).pack(side='right', padx=4)
            canvas = FigureCanvasTkAgg(fig, master=win)
            canvas.draw()
            canvas.get_tk_widget().pack(fill='both', expand=True)
            lineplot_state.update(
                {'win': win, 'axes': axes, 'fig': fig, 'canvas': canvas, 'mode': mode})
        return lineplot_state['axes'], lineplot_state['fig'], lineplot_state['canvas']

    def doColPlot(col):
        axes, fig, canvas = openOrReuseLineplot('col')
        rows = np.arange(ny)
        colors = []
        for idef, ax in zip(image_defs, axes):
            dec = idef['dec']
            if dec.ndim == 2:
                line, = ax.plot(rows, dec[:, col], label=f'col {col}')
                colors.append(line.get_color())
            else:
                plotted = [ax.plot(rows, dec[:, col, i],
                                   label=f'col {col} {ch}')[0]
                           for i, ch in enumerate(('R', 'G', 'B')[:dec.shape[2]])]
                colors.append(plotted[0].get_color())
            ax.legend(fontsize=7)
        fig.tight_layout()
        canvas.draw_idle()
        status_var.set(f'  plotted col {col}')
        return colors

    def doRowPlot(row):
        axes, fig, canvas = openOrReuseLineplot('row')
        cols_arr = np.arange(nx)
        colors = []
        for idef, ax in zip(image_defs, axes):
            dec = idef['dec']
            if dec.ndim == 2:
                line, = ax.plot(cols_arr, dec[row, :], label=f'row {row}')
                colors.append(line.get_color())
            else:
                plotted = [ax.plot(cols_arr, dec[row, :, i],
                                   label=f'row {row} {ch}')[0]
                           for i, ch in enumerate(('R', 'G', 'B')[:dec.shape[2]])]
                colors.append(plotted[0].get_color())
            ax.legend(fontsize=7)
        fig.tight_layout()
        canvas.draw_idle()
        status_var.set(f'  plotted row {row}')
        return colors

    def report_pick(col, row):
        if 0 <= row < ny and 0 <= col < nx:
            val_parts = []
            for i, idef in enumerate(image_defs):
                dec = idef['dec']
                if dec.ndim == 2:
                    val_parts.append(f'val{i+1}={dec[row, col]:.6g}')
                else:
                    val_parts.append(f'val{i+1}='
                                     + '/'.join(f'{v:.4g}' for v in dec[row, col]))
            status_var.set(f'col={col}\nrow={row}\n' + '   '.join(val_parts))
        else:
            status_var.set(f'col={col}\nrow={row}\n(out of bounds)')

    next_plot_y = [None]

    def nextPlotGeometry(win_w, win_h):
        root.update_idletasks()
        px = root.winfo_x()
        py = root.winfo_y()
        pw = root.winfo_width()
        x = max(0, min(px + pw + 10, sw - win_w))
        if next_plot_y[0] is None:
            next_plot_y[0] = py
        y = max(0, min(next_plot_y[0], sh - win_h))
        next_plot_y[0] += win_h + 10
        if next_plot_y[0] + win_h > sh:
            next_plot_y[0] = py
        return x, y

    # ---- image area: N canvases side by side, synchronized scrolling ----
    palette.update_idletasks()
    PAL_W = palette.winfo_reqwidth()
    PAL_X, PAL_Y = 10, 10
    IMG_X = PAL_X + PAL_W + 5
    IMG_Y = PAL_Y

    PLOT_WIN_W = 800
    cbar_w_total = sum(CBAR_PX if not d['is_rgb'] else 0 for d in image_defs)
    cbar_per_img = max(CBAR_PX if not d['is_rgb'] else 0 for d in image_defs)
    img_area_w = sw - IMG_X - PLOT_WIN_W - 20
    usable_h = sh - IMG_Y - 60

    # Compute viewport dimensions for both stacking orientations, pick larger area
    vw_h = min(nx, max(50, (img_area_w - n_imgs * SCROLLBAR_W - cbar_w_total) // n_imgs))
    vh_h = min(ny, max(50, usable_h))
    vw_v = min(nx, max(50, img_area_w - SCROLLBAR_W - cbar_per_img))
    vh_v = min(ny, max(50, (usable_h - n_imgs * SCROLLBAR_W) // n_imgs))
    stack_horiz = (n_imgs == 1) or (vw_h * vh_h >= vw_v * vh_v)
    viewport_w = vw_h if stack_horiz else vw_v
    viewport_h = vh_h if stack_horiz else vh_v

    all_canvases = []

    def sync_xview(*args):
        for c in all_canvases:
            c.xview(*args)

    def sync_yview(*args):
        for c in all_canvases:
            c.yview(*args)

    outer = ttk.Frame(root)
    outer.pack(fill='both', expand=True)

    for i, idef in enumerate(image_defs):
        photo = decToPhoto(idef['dec'], idef['cmap'], idef['vmin'], idef['vmax'],
                           idef['is_rgb'])
        img_frame = ttk.Frame(outer)
        img_frame.pack(side='left' if stack_horiz else 'top', fill='both', expand=True)

        if not idef['is_rgb']:
            cbar_fig = makeColorbarFig(idef['cmap'], idef['vmin'], idef['vmax'], viewport_h)
            cbar_cv  = FigureCanvasTkAgg(cbar_fig, master=img_frame)
            cbar_cv.draw()
            cbar_cv.get_tk_widget().pack(side='right', fill='y')

        sf = ttk.Frame(img_frame)
        sf.pack(side='left', fill='both', expand=True)

        ttk.Label(sf, text=f'{i+1}) {os.path.basename(idef["title"])}',
                  anchor='center').pack(side='top', fill='x', pady=(0, 1))

        xs = ttk.Scrollbar(sf, orient='horizontal')
        ys = ttk.Scrollbar(sf, orient='vertical')
        xs.pack(side='bottom', fill='x')
        ys.pack(side='right',  fill='y')

        c = tk.Canvas(sf, width=viewport_w, height=viewport_h,
                      xscrollcommand=xs.set, yscrollcommand=ys.set)
        c.pack(side='left', fill='both', expand=True)
        xs.config(command=sync_xview)
        ys.config(command=sync_yview)

        c.create_image(0, 0, anchor='nw', image=photo)
        c.image = photo
        c.config(scrollregion=(0, 0, nx, ny))
        all_canvases.append(c)

    def _wy(event):
        for c in all_canvases: c.yview_scroll(int(-1 * (event.delta / 120)), 'units')
    def _wx(event):
        for c in all_canvases: c.xview_scroll(int(-1 * (event.delta / 120)), 'units')
    def _b4(event):
        for c in all_canvases: c.yview_scroll(-1, 'units')
    def _b5(event):
        for c in all_canvases: c.yview_scroll(1, 'units')
    def _sb4(event):
        for c in all_canvases: c.xview_scroll(-1, 'units')
    def _sb5(event):
        for c in all_canvases: c.xview_scroll(1, 'units')
    for c in all_canvases:
        c.bind('<MouseWheel>',       _wy)
        c.bind('<Shift-MouseWheel>', _wx)
        c.bind('<Button-4>',         _b4)
        c.bind('<Button-5>',         _b5)
        c.bind('<Shift-Button-4>',   _sb4)
        c.bind('<Shift-Button-5>',   _sb5)

    # ---- overlay helpers (items stored as (canvas, item_id) pairs) ----
    profile_overlay_items = []
    colrow_overlay_items  = []

    def clear_overlay():
        for canvas, item in profile_overlay_items:
            canvas.delete(item)
        profile_overlay_items.clear()

    def _add_canvas_item(canvas, item, lst):
        if not lines_visible[0]:
            canvas.itemconfigure(item, state='hidden')
        lst.append((canvas, item))

    def draw_marker(col, row, color='yellow'):
        r = 5
        for canvas in all_canvases:
            item = canvas.create_oval(col - r, row - r, col + r, row + r,
                                      outline=color, width=2)
            _add_canvas_item(canvas, item, profile_overlay_items)

    def draw_profile_line(c0, r0, c1, r1, color='yellow'):
        for canvas in all_canvases:
            item = canvas.create_line(c0, r0, c1, r1, fill=color, width=1, dash=(4, 2))
            _add_canvas_item(canvas, item, profile_overlay_items)

    def draw_col_line(col, colors):
        for canvas, color in zip(all_canvases, colors):
            item = canvas.create_line(col, 0, col, ny - 1, fill=color, width=1)
            _add_canvas_item(canvas, item, colrow_overlay_items)

    def draw_row_line(row, colors):
        for canvas, color in zip(all_canvases, colors):
            item = canvas.create_line(0, row, nx - 1, row, fill=color, width=1)
            _add_canvas_item(canvas, item, colrow_overlay_items)

    def canvas_set_visible(v):
        state = 'normal' if v else 'hidden'
        for canvas, item in profile_overlay_items + colrow_overlay_items:
            canvas.itemconfigure(item, state=state)
    overlay_set_visible[0] = canvas_set_visible

    # ---- click handler (bound to all canvases) ----
    def on_canvas_click(event):
        canvas = event.widget
        col = int(canvas.canvasx(event.x))
        row = int(canvas.canvasy(event.y))

        if pick_active[0]:
            report_pick(col, row)

        elif profile_active[0]:
            if len(profile_pts) == 0:
                clear_overlay()
                profile_pts.append((row, col))
                draw_marker(col, row)
                status_var.set(f'  Profile: first point col={col} row={row}'
                               f'  — click second point')
            else:
                profile_pts.append((row, col))
                draw_marker(col, row)
                p0, p1 = profile_pts
                draw_profile_line(p0[1], p0[0], p1[1], p1[0])
                dv        = [extractProfile(d['dec'], p0[0], p0[1], p1[0], p1[1])
                             for d in image_defs]
                dist      = dv[0][0]
                vals_list = [x[1] for x in dv]
                titles    = [d['title'] for d in image_defs]
                prof_h    = 200 + 200 * n_imgs
                openProfileWindow(dist, vals_list, p0, p1, titles=titles,
                                  pos=nextPlotGeometry(800, prof_h))
                status_var.set(f'  Profile: ({p0[1]},{p0[0]}) → ({p1[1]},{p1[0]})'
                               f'  {dist[-1]:.1f} px — click to start new profile')
                profile_pts.clear()

        elif col_active[0]:
            colors = doColPlot(col)
            draw_col_line(col, colors)

        elif row_active[0]:
            colors = doRowPlot(row)
            draw_row_line(row, colors)

    for c in all_canvases:
        c.bind('<Button-1>', on_canvas_click)

    # ---- position palette at left, image window to the right ----
    if stack_horiz:
        win_w = n_imgs * (viewport_w + SCROLLBAR_W) + cbar_w_total
        win_h = viewport_h + SCROLLBAR_W
    else:
        win_w = viewport_w + SCROLLBAR_W + cbar_per_img
        win_h = n_imgs * (viewport_h + SCROLLBAR_W)

    status_lbl.config(wraplength=max(60, PAL_W - 12))
    palette.geometry(f'{PAL_W}x{win_h}+{PAL_X}+{PAL_Y}')
    root.geometry(f'{win_w}x{win_h}+{IMG_X}+{IMG_Y}')

    root.mainloop()


def main():
    parser = argparse.ArgumentParser(
        description='Display 1–3 same-size VRT or GeoTIFF images side by side.',
        epilog='Part of the utilities package.')
    parser.add_argument('files', metavar='FILE', nargs='+',
                        help='Input image(s) (.vrt, .tif, .tiff) — up to 3')
    parser.add_argument('--cmap', default='gray',
                        help='Colormap for single-band images (default: gray)')
    parser.add_argument('--vmin', type=float, default=None,
                        help='Lower clip value (default: 2nd percentile)')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Upper clip value (default: 98th percentile)')
    parser.add_argument('-vmin', type=float, default=None, dest='vmin',
                        help=argparse.SUPPRESS)
    parser.add_argument('-vmax', type=float, default=None, dest='vmax',
                        help=argparse.SUPPRESS)
    parser.add_argument('--decFactor', type=int, default=None,
                        help='Decimation factor (default: auto-fit to screen)')
    parser.add_argument('-decFactor', type=int, default=None, dest='decFactor',
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if len(args.files) > 3:
        sys.exit('showimage: at most 3 files can be displayed simultaneously')

    try:
        from osgeo import gdal
    except ImportError:
        sys.exit('osgeo.gdal not available — install gdal')

    gdal.UseExceptions()

    datasets = []
    for f in args.files:
        try:
            datasets.append(gdal.Open(f))
        except Exception as e:
            sys.exit(f'Cannot open {f}: {e}')

    sizes = [(ds.RasterXSize, ds.RasterYSize) for ds in datasets]
    if len(set(sizes)) > 1:
        msgs = [f'  {f}: {w}×{h}' for f, (w, h) in zip(args.files, sizes)]
        sys.exit('showimage: all images must have the same dimensions:\n' + '\n'.join(msgs))

    nx, ny = sizes[0]
    sw, sh = getScreenSize()

    if args.decFactor is not None:
        factor = max(1, args.decFactor)
    else:
        factor = max(1, math.ceil(nx / sw), math.ceil(ny / sh))

    image_defs = []
    for ds, f in zip(datasets, args.files):
        nb = ds.RasterCount
        print(f'{f}: {nx}×{ny} px, {nb} band(s), decimation ×{factor}')

        if nb == 1:
            dec = blockAverage(readBand(ds, 1), factor)
            vmin = args.vmin if args.vmin is not None else np.nanpercentile(dec, 2)
            vmax = args.vmax if args.vmax is not None else np.nanpercentile(dec, 98)
            is_rgb = False
        else:
            n_read = min(nb, 3)
            bands = np.stack([readBand(ds, i) for i in range(1, n_read + 1)], axis=-1)
            dec = blockAverage(bands, factor)
            if n_read == 3:
                for i in range(3):
                    lo = np.nanpercentile(dec[:, :, i], 2)
                    hi = np.nanpercentile(dec[:, :, i], 98)
                    dec[:, :, i] = np.clip(
                        (dec[:, :, i] - lo) / max(hi - lo, 1e-10), 0, 1)
                dec = np.nan_to_num(dec, nan=0.0)
                vmin = vmax = None
                is_rgb = True
            else:
                dec = dec[:, :, 0]
                vmin = args.vmin if args.vmin is not None else np.nanpercentile(dec, 2)
                vmax = args.vmax if args.vmax is not None else np.nanpercentile(dec, 98)
                is_rgb = False

        image_defs.append({
            'dec': dec,
            'title': f,
            'cmap': args.cmap,
            'vmin': vmin,
            'vmax': vmax,
            'is_rgb': is_rgb,
        })

    showImage(image_defs, sw, sh)


if __name__ == '__main__':
    main()
