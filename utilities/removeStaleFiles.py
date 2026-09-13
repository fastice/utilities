# removeStaleFiles.py
import glob
import os


def removeStaleFiles(pattern, cutoffTime, verbose=False):
    """ Remove files matching a glob pattern whose mtime predates cutoffTime.

    For clearing products left by an earlier run once a rerun has written its
    own. Everything the rerun produced is newer than cutoffTime (the time the
    run started), so whatever is still older is stale -- including output
    written in a different format, which a rerun does not overwrite: a binary
    mosaicOffsets.vx sits happily beside the mosaicOffsets.vx.tif of a later
    GeoTIFF run, and nothing but its age distinguishes it.

    Only call this after confirming the run succeeded. On a failed run nothing
    new was written, every file predates cutoffTime, and this would delete the
    products the failed run was meant to replace.

    Passing time.time() as cutoffTime removes every current match, which is
    how a caller drops bands it never wants to keep.

    Directories are left alone. Returns the list of files removed.
    """
    removed = []
    for myFile in sorted(glob.glob(pattern)):
        if not os.path.isfile(myFile):
            continue
        try:
            if os.path.getmtime(myFile) < cutoffTime:
                os.remove(myFile)
                removed.append(myFile)
        except OSError as e:
            print(f'removeStaleFiles: could not remove {myFile} ({e})')
    if verbose and removed:
        print(f'removeStaleFiles: removed {len(removed)} stale file(s) '
              f'matching {pattern}')
    return removed
