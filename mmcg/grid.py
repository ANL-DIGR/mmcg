""" Module that takes a radar object, processed using the package CMAC 2.0
a maps """

import numpy as np
import pyart

from .config import get_grid_values

def mmcg(radar, grid_shape, grid_limits, config=None,
         verbose=True, **kwargs):
    """
    Mapped Moments to a Cartesian Grid

    Parameters
    ----------
    radar : Radar
        Radar object, processed by CMAC 2.0, to be mapped to a Cartesian
        grid.
    grid_shape : 3-tuple of floats
        Number of points in the grid (z, y, x).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        z, y, x coordinates.

    Other Parameters
    ----------------
    config : str
        A string pointing to dictionaries containing values for gridding.
        These dictionaries tend to have value for grid_shape and grid_limits,
        but also all parameters found in map_gates_to_grid. This allows for
        inputing values that the user found to be best for interpolation
        and artifact removal, and call these values again from a named
        configuration.
    verbose : bool
        If True, this will display more statistics.
    kwargs : **kwargs
        Parameters found in map_gates_to_grid. For more detail:
        <https://github.com/ARM-DOE/pyart/blob/master/pyart/map/gates_to_
        grid.py#L30-L157>

    Returns
    -------
    grid : Grid
        Radar object with new CMAC added fields.

    """
    # Retrieve values from the configuration file.
    if config is not None:
        grid_config = get_grid_values(config)

        # Define gridding parameters based on configuration.
        grid_origin = grid_config['grid_origin']
        grid_origin_alt = grid_config['grid_origin_alt']
        grid_projection = grid_config['grid_projection']
        fields = grid_config['fields']
        map_roi = grid_config['map_roi']
        weighting_function = grid_config['weighting_function']
        toa = grid_config['toa']
        roi_func = grid_config['roi_func']
        constant_roi = grid_config['constant_roi']
        z_factor = grid_config['z_factor']
        xy_factor = grid_config['xy_factor']
        min_radius = grid_config['min_radius']
        h_factor = grid_config['h_factor']
        nb = grid_config['nb']
        bsp = grid_config['bsp']

        grid = pyart.map.grid_from_radars(
            radar, grid_shape, grid_limits, fields, map_roi,
            weighting_function, toa, roi_func, constant_roi, z_factor,
            xy_factor, min_radius, h_factor, nb, bsp)

    else:
        grid = pyart.map.grid_from_radars(
            radar, grid_shape, grid_limits, **kwargs)

    return grid
