# Convert triangulation polygon values into gridded netcdf

import numpy as np
import matplotlib.path


def tri2grid(grid_inc, range_strain, triangle_vertices, rot, exx, exy, eyy):
    """
    Bring delaunay 1-D quantities into the same 2-D form as the other methods

    :param grid_inc: list
    :param range_strain: list
    :param triangle_vertices: list
    :param rot: 1D array
    :param exx: 1D array
    :param exy: 1D array
    :param eyy: 1D array
    """
    lons, lats, grid = make_grid(range_strain, grid_inc);
    print("Producing gridded dataset of: Exx")
    exx_grd = find_in_triangles(triangle_vertices, exx, lons, lats, grid);
    print("Producing gridded dataset of: Exy")
    exy_grd = find_in_triangles(triangle_vertices, exy, lons, lats, grid);
    print("Producing gridded dataset of: Eyy")
    eyy_grd = find_in_triangles(triangle_vertices, eyy, lons, lats, grid);
    print("Producing gridded dataset of: Rot")
    rot_grd = find_in_triangles(triangle_vertices, rot, lons, lats, grid);
    return lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd;


def make_grid(coordbox, inc):
    """
    :param coordbox: [float, float, float, float] corresponding to [W, E, S, N]
    :type coordbox: list
    :param inc: [float, float] corresponding to [xinc, yinc]
    :type inc: list
    :returns: 1d array of lons, 1d array of lats, 2d array of zeros
    """
    lonmin = coordbox[0]
    lonmax = coordbox[1]
    latmin = coordbox[2]
    latmax = coordbox[3]
    lons = np.arange(lonmin, lonmax+0.00001, inc[0])
    lats = np.arange(latmin, latmax+0.00001, inc[1])
    grid = np.zeros((len(lats), len(lons)));
    return lons, lats, grid


def find_in_triangles(triangles, values, lons, lats, grid):
    """
    search triangle vertices for each gridpoint, then assigns that triangle's value to gridpoint

    :param triangles: list
    :param values: list
    :param lons: list
    :param lats: list
    :param grid: 2D array
    """
    val_arr = np.nan*np.ones(np.shape(grid));
    for j in range(np.shape(grid)[0]):
        for k in range(np.shape(grid)[1]):
            for i in range(len(triangles)):
                tri_lons = [triangles[i][0][0], triangles[i][1][0], triangles[i][2][0]];
                tri_lats = [triangles[i][0][1], triangles[i][1][1], triangles[i][2][1]];
                if np.min(tri_lons) < lons[k] < np.max(tri_lons):  # an extra check to speed up the process
                    if np.min(tri_lats) < lats[j] < np.max(tri_lats):
                        tripath = matplotlib.path.Path(triangles[i])
                        if tripath.contains_point((lons[k], lats[j])):
                            val_arr[j][k] = values[i];
                        break;
    return val_arr
