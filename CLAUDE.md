# CLAUDE.md — utilities

Foundational helper package used by nearly every other GrIMP Python package (`sarfunc`, `mosaicfunc`, `mosaicworkflow`, `nisargrimpworkflow`, `insarworkflow`, ...). Conventionally imported as `import utilities as u`. See the [packages CLAUDE.md](../CLAUDE.md) for pipeline context and the binary/`.dat`/VRT format contract.

## Module Map

| Module | Exports | Purpose |
|---|---|---|
| `geodat.py` | `geodat` | Polar-stereographic image geometry (origin/spacing/size) + ll↔xy transforms |
| `geoimage.py` | `geoimage` | Read/write GrIMP velocity/error/scalar binary images + GeoTIFFs |
| `geodatrxa.py` | `geodatrxa` | SAR range/azimuth geometry (`geodatNxM.in`/`.geojson`), state vectors, geocoding |
| `offsets.py` | `offsets` | Range/azimuth offset maps (`.da`/`.dr`), sigma, mask, VRT I/O |
| `lsdat.py` | `lsdat` | Landsat match `.dat` metadata + geocoding |
| `lsfit.py` | `lsfit` | Parse `lsfit` (Landsat tie-point fit) output files |
| `readImage.py` / `writeImage.py` | `readImage`, `writeImage` | Raw binary flat-file I/O |
| `readwriteLLtoRA.py` | `readLLtoRA`, `writeLLtoRAformat` | Binary I/O for the `lltora` C tool |
| `runMyThreads.py` | `runMyThreads` | Thread-pool runner with progress display |
| `pushdpopd.py` | `pushd`, `popd` | csh-style directory stack |
| `myerror.py` / `mywarning.py` / `myalert.py` | `myerror`, `mywarning`, `myalert` | Colored error/warning/info messages |
| `myPrompt.py` | `myPrompt` | y/n interactive prompt |
| `callMyProg.py` | `callMyProg` | Subprocess wrapper (stdin input list, optional logger) |
| `logger.py` | `logger` | Timestamped run log with arg dumping |
| `getWKT_PROJ.py` | `getWKT_PROJ` | EPSG → WKT1_GDAL string (via pyproj) |
| `makeMaskFromShape.py` | `makeMaskFromShape` | Rasterize a shapefile into a 0/1 mask for a `geodat` |
| `shpplot.py` | `shpplot` | Time-series / profile plotting helper (matplotlib) |
| `hsvVelCmap.py` | `hsvVelCmap` | HSV-based velocity colormap (for browse images) |
| `processProfile.py` | `processProfile` | Resample a polyline at near-uniform spacing |
| `runvel.py` | `runvel` | Thin wrapper that calls `mosaic3d -center ...` |
| `strip.py` | `strip` | Parse `#Ncol .. & ` delimited numeric data blocks |
| `dols.py` | `dols` | Run an `ls`-style shell command via csh, return token list |

## `geodat` — polar-stereographic image geometry

```python
g = u.geodat(domain='greenland')   # or 'antarctica'; sets EPSG:3413 / EPSG:3031
g.readGeodat('range.offsets.geodat')   # legacy text format
g.readGeodatFromTiff('out.vx.tif')     # derive geodat from a GeoTIFF geotransform
nx, ny = g.sizeInPixels()
x0, y0 = g.originInKm()
dx, dy = g.pixSizeInM()
x, y = g.lltoxykm(lat, lon)            # WGS84 -> PS km
xi, yi = g.xymtoImage(x, y)            # PS meters -> image pixel coords
```

Constructor accepts `wkt=` (WKT or `EPSG:nnnn` string) to override the domain-derived projection. `writeGeodat(path)` writes the legacy 3-line `# 2 / nx ny / dx dy / x0 y0 / &` format.

## `geoimage` — velocity/error/scalar image I/O

Type-driven via `geoType` in `{'scalar', 'velocity', 'velocityRA', 'error'}`. `_TYPE_CONFIG` (module-level dict) maps each type to component attribute names, file suffixes, and the geodat-sidecar suffix:

| geoType | components | suffixes | geodat suffix | magnitude attr |
|---|---|---|---|---|
| `scalar` | `x` | `''` | `.geodat` | none |
| `velocity` | `vx`, `vy` | `.vx`, `.vy` | `.vx.geodat` | `v` |
| `velocityRA` | `vr`, `va` | `.vr`, `.va` | `.vr.geodat` | `v` |
| `error` | `ex`, `ey` | `.ex`, `.ey` | `.ex.geodat` | `e` |

```python
img = u.geoimage(verbose=False)
img.readData('Vel-2022-01-01.2022-01-31/release/myfile', geoType='velocity',
              epsg=3413)             # reads myfile.vx/.vy + myfile.vx.geodat
img.xyGrid()                          # populate img.xGrid/img.yGrid (km)
vx, vy, v = img.interpGeo(xq, yq)     # bilinear interp at query points (km)
img.writeMyTiff('out/myfile', epsg=3413, driverName='COG')  # writes .vx.tif, .vy.tif, .v.tif
img.writeData('out/myfile', dType='>f4')  # writes binary flat files + .geodat sidecars
```

`readFiles` replaces values `<= -2e9` (the GrIMP no-data sentinel) with `NaN` on read; `writeData`/`writeCloudOptGeo` reverse this (`NaN -> -2e9` or the type-specific `_NO_DATA` value) on write. `noV=True` in `writeMyTiff` skips writing the magnitude (`.v`/`.e`) band.

## `offsets` — range/azimuth offset maps

Models `azimuth.offsets`/`range.offsets` (or `*.da`/`*.dr`) plus their `.dat`, `.sa`/`.sr` sigma, `.mt` matchType, and `.mask` siblings. File-name derivation rules live in `offsetFileNames`, `datFileName`, `sigmaFileName`, `matchTypeFileName`, `maskFileName` — each tries a default based on `fileRoot` and accepts an explicit override.

```python
off = u.offsets(fileRoot='offsets/azimuth.offsets', verbose=False)
rg, az = off.readOffsets()            # -> float arrays, -2e9 = no data
off.readSigma()
off.readMask()
lat, lon = off.getLatLon()
heading = off.computeHeading()        # needs geodatrxaFile + lat/lon
off.writeOffsets(byteOrder='MSB')     # writes .dr/.da + .dat sidecars
```

If `vrtFile` is set (or a `.vrt` exists alongside the `.dat`), all reads route through `readVrt()` instead of raw binary — this is the path used for NISAR products where a combined `offsets.range-azimuth.vrt` carries both data and metadata (`a0`, `r0`, `deltaA`, `deltaR`, `sigmaRange`/`sigmaStreaks` -> `rgErr`/`azErr`, `ByteOrder`). `-1.99e9` is the validity threshold used by `areValid()`/`notValid()`.

## `geodatrxa` — SAR range/azimuth geometry

Reads either legacy `geodatNxM.in` or (preferred) `geodatNxM.geojson`. `readFile()` auto-substitutes `.geojson` for `.in` if present unless `forceIn=True`. Carries orbit state vectors (`position`/`velocity`/`stateTime`), corner coordinates, look direction/pass type, and provides:

- `isSouth()`, `isRightLooking()`, `isDescending()`/`isAscending()`
- `nearRangem()`/`farRangem()`/`centerRangem()`/`satelliteAltm()`
- `interpPos(t)`/`interpVel(t)` — cubic-interpolated state vectors
- `llzPtToRA(lat, lon, z)` — iterative geocoding to range/azimuth/time
- `thetaCrad()`/`thetaCActualrad()` — incidence angle from geometry
- `writeGeodatFile(path)` — writes legacy `.in` format from a geojson-loaded object

## `runMyThreads` and `pushd`/`popd`

`runMyThreads(threads, maxThreads, message, delay=0.2, prompt=False, quiet=False)` is the standard parallelism primitive: takes a list of unstarted `threading.Thread` objects, runs up to `maxThreads` concurrently, prints a live progress line, and always `os.chdir(home)` before starting each thread (so threads that `pushd`/`popd` don't corrupt the working directory for siblings). Calls `myalert('Threads Done')` on completion. Used by `mosaicworkflow.setupquarters` for per-sector `mosaic3d` runs.

`pushd(directory)` / `popd()` maintain a module-level directory stack (not thread-safe — `runMyThreads` resets cwd between thread starts to mitigate this).

## Error/warning helpers

```python
u.myerror(message)     # red banner, sys.exit()
u.mywarning(message)   # yellow banner, continue
u.myalert(message)     # cyan banner, continue (status messages)
```

`myerror` accepts `myLogger=` to also call `logger.logError()` before exiting.

## Notes

- **`-2e9` sentinel**: the GrIMP no-data value for float32 images. `geoimage.readFiles` converts `<= -2e9` to `NaN` on read and back on write. `offsets` uses `-1.99e9` as its "valid" threshold (`areValid`/`notValid`) — slightly looser to tolerate roundoff.
- **Byte order**: `readImage`/`writeImage` accept `'>f4'` etc.; the `'>'` prefix triggers an explicit byteswap to native order on read (and back to MSB on write). This matches the unconditional byteswap in GIT64's `freadBS`/`fwriteBS`.
- **geodat vs geodatrxa**: `geodat` is for *output mosaic* geometry (polar stereographic, km-based). `geodatrxa` is for *input SAR product* range/azimuth geometry (single-look pixel size, orbit state vectors). Don't confuse the two — `offsets.geodatrxa` (an attribute) holds a `geodatrxa` instance for the offset's source SAR geometry.
- **`.in` vs `.geojson`**: `geodatrxa.readFile()` prefers `.geojson` if both exist for the same base name; `offsets.geodatrxaFile` and `sarDB`/`sensorDefinitions.geodatName()` (in `sarfunc`) follow the same `geodatType` convention (`'in'` or `'geojson'`).
- `shpplot` requires `PIL`/matplotlib at import time for some submodules (`makeMaskFromShape` imports `PIL.Image`/`ImageDraw`) — wrapped in `try/except` with a warning print if unavailable.
- `dols()` shells out via `/bin/csh` explicitly — relevant if invoked from a non-csh environment.
- Several modules (`runvel.py`, `myPrompt.py`, `callMyProg.py`) are thin/legacy wrappers from early in the project; `runvel` hardcodes a `mosaic3d -center` invocation and is largely superseded by `mosaicworkflow.setupquarters`'s own threaded calls.
