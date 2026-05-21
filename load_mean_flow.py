# -*- coding: utf-8 -*-

import numpy as np

from scipy.io import loadmat
from scipy.spatial import Delaunay

from interpolation import (
    build_linear_interpolation_matrix,
    build_nearest_interpolation_matrix,
    apply_interpolation
)


def load_mean_flow(ni, nj, ni_ansys, nj_ansys):
    

    """
    Loads the structured Pointwise mesh and the mean flow solution,
    then constructs the interpolation operators required to map
    the ANSYS flow fields onto the structured mesh.

    The interpolated mean fields

        p_mean, u_mean, v_mean

    are used to compute the fluctuating flow quantities.
    """

    with open('aoa_18.x', 'r') as f:

        grid = np.fromstring(
            f.read(),
            sep=' '
        )

    print("Number of elements:", grid.size)

    x_coords_flat = grid[:grid.size // 2]
    y_coords_flat = grid[grid.size // 2:]

    expected_points = ni * nj

    if x_coords_flat.size > expected_points:

        extra = x_coords_flat.size - expected_points

        print('Removing extra points = ', extra)

        x_coords_flat = x_coords_flat[:-extra]
        y_coords_flat = y_coords_flat[:-extra]

    x_grid = np.reshape(
        x_coords_flat,
        (ni, nj),
        order='F'
    )

    y_grid = np.reshape(
        y_coords_flat,
        (ni, nj),
        order='F'
    )

    struct = loadmat('mean_puv_aoa18.mat')

    data = np.reshape(

        struct['mean_puv_aoa'],
        (ni_ansys, nj_ansys, 6),
        order='F'

    )

    x_source = data[:, :, 1]
    y_source = data[:, :, 2]

    p_unordered = data[:, :, 3]
    u_unordered = data[:, :, 4]
    v_unordered = data[:, :, 5]

    source_points = np.column_stack((

        x_source.flatten(order='F'),
        y_source.flatten(order='F')

    ))

    tri = Delaunay(source_points)

    W_linear = build_linear_interpolation_matrix(

        tri,
        source_points,
        x_grid,
        y_grid

    )

    W_nearest = build_nearest_interpolation_matrix(

        source_points,
        x_grid,
        y_grid

    )

    p_mean = apply_interpolation(
        W_linear,
        W_nearest,
        p_unordered,
        ni,
        nj
    )

    u_mean = apply_interpolation(
        W_linear,
        W_nearest,
        u_unordered,
        ni,
        nj
    )

    v_mean = apply_interpolation(
        W_linear,
        W_nearest,
        v_unordered,
        ni,
        nj
    )

    return (

        W_linear,
        W_nearest,

        p_unordered,
        u_unordered,
        v_unordered,

        p_mean,
        u_mean,
        v_mean

    )