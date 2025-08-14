import dask.array.core
import xarray
import matplotlib
import colorcet as cc
import numpy as np

def best_datatree_scale(dt: xarray.DataTree, max_dim: int) -> str:
    '''Returns the scale of the DataTree level whose x and y dimensions are 
    less than or equal to max_dim. If all scales have larger dimensions, 
    returns the lowest resolution scale.'''
    assert isinstance(dt, xarray.DataTree), f"Expected xarray.DataTree, got {type(dt)}"
    for scale in dt:  # check successively lower resolutions (scale0, scale1, ...)
        if dt[scale].sizes['x'] <= max_dim and dt[scale].sizes['y'] <= max_dim: 
            return scale
    # All scales are too large, so return the lowest resolution one
    return list(dt.keys())[-1]

def rel_fontsize(fig, frac_height):
    height_in = fig.get_size_inches()[1]      # Figure height in inches
    return int(72 * frac_height * height_in)  # Matplotlib's fontsize is given in points (1pt = 1/72 inch)

def glasbey_cmap(n):
    base = cc.cm['glasbey']
    colors = base(np.linspace(0, 1, n))
    return matplotlib.colors.ListedColormap(colors)