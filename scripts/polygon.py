"""
Replacement for descartes PolygonPatch.

See https://github.com/geopandas/geopandas/issues/1039#issuecomment-748625852

Source: https://github.com/geopandas/geopandas/blob/831108dcb369517e311cd0ee3821db1f2fb7122d/geopandas/plotting.py#L100-L123
"""

import numpy as np

def PolygonCollection(polygons, **kwargs):
    """
    Constructs a matplotlib PatchCollection from polygon geometries.
    """
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    from matplotlib.collections import PatchCollection
    
    return PatchCollection([
        PathPatch(Path.make_compound_path(
            Path(np.asarray(polygon.exterior.coords)[:, :2]),
            *[ Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors ]))
        for polygon in polygons.geoms if not polygon.is_empty ],
                           **kwargs)


