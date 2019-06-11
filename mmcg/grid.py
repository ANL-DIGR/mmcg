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
        Whether or not to map fields in origin dBZ units in linear or
        logarithmic units. Default is True, fields are mapped in linear units
        and then converted to logarithmic units after mapping has taken place.
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
        total_pow = 10.0**(radar.fields['total_power']['data']/10.0)
        corr_z = 10.0**(radar.fields['corrected_reflectivity']['data']/10.0)
        radar.add_field_like(
            'reflectivity', 'reflectivity', z_lin, replace_existing=True)
        radar.add_field_like(
            'total_power', 'total_power', total_pow, replace_existing=True)
        radar.add_field_like(
            'corrected_reflectivity', 'corrected_reflectivity',
            corr_z, replace_existing=True)

    # Retrieve values from the configuration file.
    if config is not None:
        grid_config = get_grid_values(config)

        grid = pyart.map.grid_from_radars(
            radar, grid_shape, grid_limits, **grid_config)
        if 'gate_id' in radar.fields.keys():
            if 'fields' in grid_config:
                grid_config.pop('fields')
            if 'weighting_function' in grid_config:
                grid_config.pop('weighting_function')
            grid_id = pyart.map.grid_from_radars(
                radar, grid_shape, grid_limits, fields=['gate_id'],
                weighting_function='NEAREST', **grid_config)
            gate_data = grid_id.fields['gate_id']['data']
            grid.fields['gate_id']['data'] = gate_data
            grid.fields['gate_id'].update({
                'comment': 'This gate id field has been mapped to a '
                           'Cartesian grid using nearest neighbor. This '
                           'may differ from the mapping method used '
                           'in the other fields'})
            del grid_id

    else:
        grid = pyart.map.grid_from_radars(
            radar, grid_shape, grid_limits, **kwargs)
        if 'gate_id' in radar.fields.keys():
            if 'fields' in kwargs:
                kwargs.pop('fields')
            if 'weighting_function' in kwargs:
                kwargs.pop('weighting_function')
            grid_id = pyart.map.grid_from_radars(
                radar, grid_shape, grid_limits, fields=['gate_id'],
                weighting_function='NEAREST', **kwargs)
            gate_data = grid_id.fields['gate_id']['data']
            grid.fields['gate_id']['data'] = gate_data
            grid.fields['gate_id'].update({
                'comment': 'This gate id field has been mapped to a '
                           'Cartesian grid using nearest neighbor. This '
                           'may differ from the mapping method used '
                           'in the other fields'})
            del grid_id

    # Convert reflectivity back into logarithmic units and add a comment
    # in the field dictionary explaining the method used.
    if z_linear_interp:
        ref_log_grid = 10.0*(np.log10(grid.fields['reflectivity']['data']))
        grid.fields['reflectivity']['data'] = ref_log_grid
        grid.fields['reflectivity'].update({
            'comment': 'This reflectivity field was interpolated linearly and '
                       'then converted to logarithmic units. Using linear '
                       'units during interpolation allows for the retention '
                       'of storm structure and gives a more realistic '
                       'estimation of convection and more.'})

        tot_log_grid = 10.0*(np.log10(grid.fields['total_power']['data']))
        grid.fields['total_power']['data'] = tot_log_grid
        grid.fields['total_power'].update({
            'comment': 'This total power field was interpolated linearly and '
                       'then converted to logarithmic units. Using linear '
                       'units during interpolation allows for the retention '
                       'of storm structure and gives a more realistic '
                       'estimation of convection and more.'})

        corr_z_log_grid = 10.0*(
            np.log10(grid.fields['corrected_reflectivity']['data']))
        grid.fields['corrected_reflectivity']['data'] = corr_z_log_grid
        grid.fields['corrected_reflectivity'].update({
            'comment': 'This corrected_reflectivity field was interpolated '
                       'linearly and then converted to logarithmic units. '
                       'Using linear units during interpolation allows for '
                       'the retention of storm structure and gives a more '
                       'realistic estimation of convection and more.'})
    return grid
