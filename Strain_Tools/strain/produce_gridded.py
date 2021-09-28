# Convert triangulation polygon values into gridded netcdf

import numpy as np
import matplotlib.path

def tri2grid(lons, lats, triangle_vertices, rot, exx, exy, eyy):
    """
    Bring delaunay 1-D quantities into the same 2-D form as the other methods

    :param lons: list
    :param lats: list
    :param triangle_vertices: list
    :param rot: 1D array
    :param exx: 1D array
    :param exy: 1D array
    :param eyy: 1D array
    """
    print("Producing gridded dataset of: Exx")
    exx_grd = find_in_triangles(triangle_vertices, exx, lons, lats);
    print("Producing gridded dataset of: Exy")
    exy_grd = find_in_triangles(triangle_vertices, exy, lons, lats);
    print("Producing gridded dataset of: Eyy")
    eyy_grd = find_in_triangles(triangle_vertices, eyy, lons, lats);
    print("Producing gridded dataset of: Rot")
    rot_grd = find_in_triangles(triangle_vertices, rot, lons, lats);
    return rot_grd, exx_grd, exy_grd, eyy_grd;


def find_in_triangles(triangles, values, lons, lats):
    """
    search triangle vertices for each gridpoint, then assigns that triangle's value to gridpoint

    :param triangles: list
    :param values: list
    :param lons: list
    :param lats: list
    """
    val_arr = np.nan*np.ones((len(lats), len(lons)));

    for i, triangle in enumerate(triangles):
        tri_lons = [triangle[0][0], triangle[1][0], triangle[2][0]];
        tri_lats = [triangle[0][1], triangle[1][1], triangle[2][1]];
        above_min_lon = np.where(lons >= np.min(tri_lons))[0];
        below_max_lon = np.where(lons <= np.max(tri_lons))[0];
        fitting_lon = list(set(above_min_lon) & set(below_max_lon));
        above_min_lat = np.where(lats >= np.min(tri_lats))[0];
        below_max_lat = np.where(lats <= np.max(tri_lats))[0];
        fitting_lat = list(set(above_min_lat) & set(below_max_lat));

        tripath = matplotlib.path.Path(triangle)
        for j in range(len(fitting_lat)):
            for k in range(len(fitting_lon)):
                if tripath.contains_point((lons[fitting_lon[k]], lats[fitting_lat[j]])):
                    x_idx = np.where(lons == lons[fitting_lon[k]])[0];
                    y_idx = np.where(lats == lats[fitting_lat[j]])[0];
                    val_arr[y_idx[0]][x_idx[0]] = values[i];

    return val_arr;
