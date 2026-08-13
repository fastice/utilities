import glob


def globSort(pattern):
    '''
    glob.glob() does not guarantee any particular order -- it returns
    whatever order the OS's directory listing (readdir) happens to report,
    which depends on filesystem internals (e.g. hashing/B-tree storage of
    directory entries), not creation time or alphabetical order. Code that
    relies on glob() result order (e.g. picking [0], or zipping two
    separately-globbed lists together by position) can silently misbehave
    in a way that varies directory-by-directory.

    Parameters
    ----------
    pattern : str
        Glob pattern, passed directly to glob.glob().

    Returns
    -------
    list
        Sorted list of matching paths.
    '''
    return sorted(glob.glob(pattern))
