import glob


def globOffsetProducts(pattern):
    '''
    Sorted data-product roots matching a raw-name glob, finding the GeoTIFF
    form as well as the raw one.

    In tiff mode a product like <root>.interp.da is written only as
    <root>.interp.da.tif, so a glob for the raw name finds nothing and callers
    that do glob(...)[0] raise IndexError, or -- worse -- silently drop the
    product. This matches <pattern>.tif too and strips the .tif again, so the
    caller keeps handling the raw-style root name it always did.

    Returning the de-.tif'd root is deliberate: callers do string surgery on
    the result (e.g. .replace('.da', '.vrt')) and write it into mosaic input
    lists, where the C side resolves it through checkForOffsetsVrt()
    (mosaicSource/common/readOffsets.c) -- which finds the .vrt whether or not
    the raw file still exists. Leaving the .tif on would produce names like
    <root>.interp.vrt.tif.

    Parameters
    ----------
    pattern : str
        Glob pattern naming the raw product, passed to glob.glob().

    Returns
    -------
    list
        Sorted list of matching product roots, without any .tif extension.
    '''
    names = set(glob.glob(pattern))
    names |= {p[:-len('.tif')] for p in glob.glob(f'{pattern}.tif')}
    return sorted(names)
