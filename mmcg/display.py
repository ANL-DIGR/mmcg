""" Code that plots fields from the CMAC radar object. """

import os
from datetime import datetime
import operator

import netCDF4
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyart

from pyart.graph.common import (
    generate_grid_name, generate_grid_time_begin)

from .config import get_plot_values

plt.switch_backend('agg')


def quicklooks(grid, config, image_directory=None,
               dd_lobes=True):
    """
    Quicklooks, images produced with regards to CMAC

    Parameter
    ---------
    grid : Grid
        Grid object that has CMAC applied to it.
    config : str
        A string of the radar name found from config.py that contains values
        for plotting, specific to the grid created from that radar.

    Optional Parameters
    -------------------
    image_directory : str
        File path to the image folder of which to save the CMAC images. If no
        image file path is given, image path defaults to users home directory.
    dd_lobes : bool
        Plot DD lobes between radars if dd_lobes is True.

    """
    if image_directory is None:
        image_directory = os.path.expanduser('~')

    grid_start_date = netCDF4.num2date(
        grid.time['data'][0], grid.time['units'])

    # Retrieve the plot parameter values based on the radar.
    plot_config = get_plot_values(config)

    save_name = plot_config['save_name']
    date_string = datetime.strftime(grid_start_date, '%Y%m%d.%H%M%S')
    combined_name = '_grid' + '.' + save_name + '.' + date_string

    min_lat = plot_config['min_lat']
    max_lat = plot_config['max_lat']
    min_lon = plot_config['min_lon']
    max_lon = plot_config['max_lon']

    # Creating a plot of reflectivity before CMAC.
    lal = np.arange(min_lat, max_lat+.2, .2)
    lol = np.arange(min_lon, max_lon+.2, .2)

    if dd_lobes:
        grid_lat = np.arange(min_lat, max_lat, 0.01)
        grid_lon = np.arange(min_lon, max_lon, 0.01)

        facility = plot_config['facility']
        if facility == 'I4':
            dms_radar1_coords = [plot_config['site_i4_dms_lon'],
                                 plot_config['site_i4_dms_lat']]
            dms_radar2_coords = [plot_config['site_i5_dms_lon'],
                                 plot_config['site_i5_dms_lat']]
        elif facility == 'I5':
            dms_radar1_coords = [plot_config['site_i5_dms_lon'],
                                 plot_config['site_i5_dms_lat']]
            dms_radar2_coords = [plot_config['site_i4_dms_lon'],
                                 plot_config['site_i4_dms_lat']]
        elif facility == 'I6':
            dms_radar1_coords = [plot_config['site_i6_dms_lon'],
                                 plot_config['site_i6_dms_lat']]
            dms_radar2_coords = [plot_config['site_i4_dms_lon'],
                                 plot_config['site_i4_dms_lat']]

        dec_radar1 = [_dms_to_decimal(
            dms_radar1_coords[0][0], dms_radar1_coords[0][1],
            dms_radar1_coords[0][2]), _dms_to_decimal(
                dms_radar1_coords[1][0], dms_radar1_coords[1][1],
                dms_radar1_coords[1][2])]
        dec_radar2 = [_dms_to_decimal(
            dms_radar2_coords[0][0], dms_radar2_coords[0][1],
            dms_radar2_coords[0][2]), _dms_to_decimal(
                dms_radar2_coords[1][0], dms_radar2_coords[1][1],
                dms_radar2_coords[1][2])]

        bca = _get_bca(dec_radar2[0], dec_radar2[1], dec_radar1[0],
                       dec_radar1[1], grid_lon, grid_lat)
        grid_lon, grid_lat = np.meshgrid(grid_lon, grid_lat)

    level = plot_config['level']

    # Plot of the raw reflectivity from the radar.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('reflectivity', level=level,
                      vmin=-8, vmax=64, mask_outside=False,
                      cmap=pyart.graph.cm_colorblind.HomeyerRainbow)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')

    plt.savefig(
        image_directory
        + '/reflectivity' + combined_name + '.png')
    plt.close()

    # Four panel plot of gate_id, velocity_texture, reflectivity, and
    # cross_correlation_ratio.
    cat_dict = {}
    print('##')
    print('## Keys for each gate id are as follows:')
    for pair_str in radar.fields['gate_id']['notes'].split(','):
        print('##   ', str(pair_str))
        cat_dict.update({pair_str.split(':')[1]:int(pair_str.split(':')[0])})
    sorted_cats = sorted(cat_dict.items(), key=operator.itemgetter(1))
    cat_colors = {'rain': 'green',
                  'multi_trip': 'red',
                  'no_scatter': 'gray',
                  'snow': 'cyan',
                  'melting': 'yellow'}
    lab_colors = ['red', 'cyan', 'grey', 'green', 'yellow']
    if 'xsapr_clutter' in radar.fields.keys():
        cat_colors['clutter'] = 'black'
        lab_colors = np.append(lab_colors, 'black')
    lab_colors = [cat_colors[kitty[0]] for kitty in sorted_cats]
    cmap = matplotlib.colors.ListedColormap(lab_colors)

    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[15, 10])
    plt.subplot(2, 2, 1)
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('gate_id', level=level, min_lon=min_lon,
                      cmap=cmap, vmin=0, vmax=5)

    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')

    cbax = plt.gca()
    if 'xsapr_clutter' in radar.fields.keys():
        tick_locs = np.linspace(
            0, len(sorted_cats) - 2, len(sorted_cats)) + 0.5
    else:
        tick_locs = np.linspace(
            0, len(sorted_cats) - 1, len(sorted_cats)) + 0.5
    display.cbs[-1].locator = matplotlib.ticker.FixedLocator(tick_locs)
    catty_list = [sorted_cats[i][0] for i in range(len(sorted_cats))]
    display.cbs[-1].formatter = matplotlib.ticker.FixedFormatter(catty_list)
    display.cbs[-1].update_ticks()
    plt.subplot(2, 2, 2)
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('reflectivity', level=level, vmin=-8, vmax=64,
                      cmap=pyart.graph.cm_colorblind.HomeyerRainbow)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')

    plt.subplot(2, 2, 3)
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('velocity_texture', level=level, vmin=0, vmax=14,
                      title=_generate_title(
                         grid, 'velocity_texture', level),
                      cmap=pyart.graph.cm.NWSRef)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca, latlon='True',
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.subplot(2, 2, 4)
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('cross_correlation_ratio', level=level, vmin=.5,
                      vmax=1, cmap=pyart.graph.cm.Carbone42)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/cmac_four_panel_plot' + combined_name + '.png')
    plt.close()

    # Creating a plot with reflectivity corrected with attenuation.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lon, lat_lines=lal)
    display.plot_grid('attenuation_corrected_reflectivity', level=level,
                      vmin=0, vmax=60.,
                      title=_generate_title(
                          grid, 'attenuation_corrected_reflectivity',
                          level),
                      cmap=pyart.graph.cm_colorblind.HomeyerRainbow)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/attenuation_corrected_reflectivity' + combined_name + '.png')
    plt.close()

    # Creating a plot of specific attenuation.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('specific_attenuation', level=level,
                      vmin=0, vmax=1.0)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/specific_attenuation' + combined_name + '.png')
    plt.close()

    # Creating a plot of corrected differential phase.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('corrected_differential_phase', level=level,
                      title=_generate_title(
                          grid, 'corrected_differential_phase',
                          level))
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/corrected_differential_phase' + combined_name + '.png')
    plt.close()

    # Creating a plot of corrected specific differential phase.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('corrected_specific_diff_phase', sweep=sweep,
                      vmin=0, vmax=6, title=_generate_title(
                          grid, 'corrected_specific_diff_phase',
                          level))
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/corrected_specific_diff_phase' + combined_name + '.png')
    plt.close()

    # Creating a plot with region dealias corrected velocity.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('corrected_velocity', level=level,
                      cmap=pyart.graph.cm.NWSVel, vmin=-30,
                      vmax=30)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/corrected_velocity' + combined_name + '.png')
    plt.close()

    # Creating a plot of rain rate A
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lones=lol, lat_lines=lal)
    display.plot_grid('rain_rate_A', level=level, vmin=0, vmax=120)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/rain_rate_A' + combined_name + '.png')
    plt.close()

    # Creating a plot of filtered corrected differential phase.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lon, lat_lines=lat)
    display.plot_grid('filtered_corrected_differential_phase', level=level,
                      title=_generate_title(
                          grid, 'filtered_corrected_differential_phase',
                          level),
                      cmap=pyart.graph.cm.Theodore16)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/filtered_corrected_differential_phase' + combined_name + '.png')
    plt.close()

    # Creating a plot of filtered corrected specific differential phase.
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[12, 8])
    display.plot_basemap(lon_lines=lol, lat_lines=lal)
    display.plot_grid('filtered_corrected_specific_diff_phase', level=level,
                      title=_generate_title(
                          grid, 'filtered_corrected_specific_diff_phase',
                          level),
                      cmap=pyart.graph.cm.Theodore16)
    if dd_lobes:
        plt.contour(grid_lon, grid_lat, bca,
                    levels=[np.pi/6, 5*np.pi/6], linewidths=2,
                    colors='k')
    plt.savefig(
        image_directory
        + '/filtered_corrected_specific_diff_phase' + combined_name + '.png')
    plt.close()


def _generate_title(grid, field, level):
    """ Generates a title for each plot. """
    time_str = generate_grid_time_begin(grid).isoformat() + 'Z'
    height = grid.z['data'][level] / 1000.
    l1 = "%s %.1f km %s " % (generate_grid_name(grid), height,
                             time_str)
    field_name = str(field)
    field_name = field_name.replace('_', ' ')
    field_name = field_name[0].upper() + field_name[1:]
    return line_one + '\n' + field_name


def _get_bca(rad1_lon, rad1_lat, rad2_lon, rad2_lat,
             grid_lon, grid_lat):
    # Beam crossing angle needs cartesian coordinate.
    p = ccrs.PlateCarree()
    p = p.as_geocentric()
    rad1 = p.transform_points(ccrs.PlateCarree().as_geodetic(),
                              np.array(rad1_lon),
                              np.array(rad1_lat))
    rad2 = p.transform_points(ccrs.PlateCarree().as_geodetic(),
                              np.array(rad2_lon),
                              np.array(rad2_lat))
    grid_lon, grid_lat = np.meshgrid(grid_lon, grid_lat)
    grid = p.transform_points(ccrs.PlateCarree().as_geodetic(),
                              grid_lon, grid_lat,
                              np.zeros(grid_lon.shape))

    # Create grid with Radar 1 in center.
    x = grid[:, :, 0] - rad1[0, 0]
    y = grid[:, :, 1] - rad1[0, 1]
    rad2 = rad2 - rad1
    a = np.sqrt(np.multiply(x, x) + np.multiply(y, y))
    b = np.sqrt(pow(x - rad2[0, 0], 2) + pow(y - rad2[0, 1], 2))
    c = np.sqrt(rad2[0, 0] * rad2[0, 0] + rad2[0, 1] * rad2[0, 1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x - rad2[0, 1]) / b)
    return np.arccos((a*a + b*b - c*c) / (2*a*b))


def _dms_to_decimal(degrees, minutes, seconds):
    if degrees > 0:
        return degrees + minutes/60 + seconds/3600
    else:
        return degrees - minutes/60 - seconds/3600
