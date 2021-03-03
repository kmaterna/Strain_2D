# Convert triangulation polygon values into gridded netcdf

import numpy as np
import matplotlib.path


def tri2grid(grid_inc, range_strain,  triangle_vertices, rot, exx, exy, eyy):
    # steps to bring delaunay 1-D quantities into the same 2-D form as the other methods
    lons, lats, grid = make_grid(range_strain, grid_inc);
    print("Producing gridded dataset of Exx")
    exx_grd = find_in_triangles(triangle_vertices, exx, lons, lats, grid);
    print("Producing gridded dataset of: Exy")
    exy_grd = find_in_triangles(triangle_vertices, exy, lons, lats, grid);
    print("Producing gridded dataset of: Eyy")
    eyy_grd = find_in_triangles(triangle_vertices, eyy, lons, lats, grid);
    print("Producing gridded dataset of: Rot")
    rot_grd = find_in_triangles(triangle_vertices, rot, lons, lats, grid);
    return lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd;


# makes grid for delaunay
def make_grid(coordbox, inc):
    # coordbox is [float, float, float, float] [W, E, S, N]
    # inc is a float
    # return value is a 2d array of zeros
    lonmin = coordbox[0]
    lonmax = coordbox[1]
    latmin = coordbox[2]
    latmax = coordbox[3]
    lons = np.arange(lonmin, lonmax+0.00001, inc[0])
    lats = np.arange(latmin, latmax+0.00001, inc[1])
    grid = np.zeros((len(lats), len(lons)));
    return lons, lats, grid


# searches path created by triangle vertices for each gridpoint, then assigns that triangle's value to the gridpoint
def find_in_triangles(triangles, values, lons, lats, grid):
    val_arr = np.nan*np.ones(np.shape(grid));
    for j in range(np.shape(grid)[0]):
        for k in range(np.shape(grid)[1]):
            for i in range(len(triangles)):
                tripath = matplotlib.path.Path(triangles[i])
                if tripath.contains_point((lons[k], lats[j])):
                    val_arr[j][k] = values[i];
                    break;
    return val_arr

