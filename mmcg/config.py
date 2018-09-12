"""
mmcg.config
===========
MMCG Configuration.

    get_grid_values
    get_plot_values

"""

from .default_config import _DEFAULT_GRID_VALUES,  _DEFAULT_PLOT_VALUES


def get_grid_values(radar):
    """
    Return the values specific to a grid for mapping the radar data to
    Cartesian coordinates.
    """
    return _DEFAULT_GRID_VALUES[radar].copy()


def get_plot_values(radar):
    """
    Return the values specific to a grid for plotting the grid fields.
    """
    return _DEFAULT_PLOT_VALUES[radar].copy()
