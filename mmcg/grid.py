""" Module that takes a radar object, processed using the package CMAC 2.0,
a maps the data to a Cartesian grid. """

import numpy as np
import pyart

from .config import get_grid_values

def mmcg(radar, grid_shape, grid_limits, z_linear_interp=True,
         config=None, **kwargs):
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
    z_linear_interp : bool
        Whether or not to map the reflectivity in linear or logarithmic
        units. Default is True, reflectivity is map in linear units and
        then converted to logarithmic units after mapping has taken place.
    config : str
        A string pointing to dictionaries containing values for gridding.
        These dictionaries tend to have value for grid_shape and grid_limits,
        but also all parameters found in map_gates_to_grid. This allows for
        inputing values that the user found to be best for interpolation
        and artifact removal, and call these values again from a named
        configuration.
    kwargs : **kwargs
        Parameters found in map_gates_to_grid. For more detail:
        <https://github.com/ARM-DOE/pyart/blob/master/pyart/map/gates_to_
        grid.py#L30-L157>

    Returns
    -------
    grid : Grid
        Radar object with new CMAC added fields.

    """
    # Replace the raw reflectivity field with a linear reflectivity field
    # that will be used in mapping.
    if z_linear_interp:
        z_lin = 10.0**(radar.fields['reflectivity']['data']/10.0)
        radar.add_field('reflectivity', z_lin, replace_existing=True)

    # Retrieve values from the configuration file.
    if config is not None:
        grid_config = get_grid_values(config)

        grid = pyart.map.grid_from_radars(
            radar, grid_shape, grid_limits, **grid_config)

    else:
        grid = pyart.map.grid_from_radars(
            radar, grid_shape, grid_limits, **kwargs)

    # Convert reflectivity back into logarithmic units and add a comment
    # in the field dictionary explaining the method used.
    if z_linear_interp:
        ref_lin_grid = 10.0*(np.log10(grid.fields['reflectivity']['data']))
        grid.fields['reflectivity']['data'] = ref_lin_grid
        grid.fields['reflectivity'].update({
            'comment': 'This reflectivity field was interpolated linearly and '
                       'then converted to logarithmic units. Using linear '
                       'units during interpolation allows for the retention '
                       'of storm structure and gives a more realistic '
                       'estimation of convection and more.'})
    return grid
