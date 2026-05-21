# -*- coding: utf-8 -*-

import numpy as np

from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix


# ============================================================
# ------------ BUILD LINEAR INTERPOLATION MATRIX -------------
# ============================================================

def build_linear_interpolation_matrix(
    tri,
    source_points,
    target_x,
    target_y
):
    
    """
   Builds the sparse linear interpolation operator W_linear
   using barycentric coordinates obtained from the
   Delaunay triangulation of the unstructured mesh.

   The interpolation is written as

       q_ordered = W_linear q_unordered

   allowing fast repeated interpolation of flow variables
   onto the structured Pointwise mesh.
   """

    target_points = np.column_stack((

        target_x.flatten(order='F'),
        target_y.flatten(order='F')

    ))

    simplices = tri.find_simplex(target_points)

    valid = simplices >= 0

    valid_indices = np.where(valid)[0]

    simplices_valid = simplices[valid]

    vertices = tri.simplices[simplices_valid]

    transform = tri.transform[simplices_valid]

    delta = (
        target_points[valid]
        - transform[:, 2]
    )

    bary = np.einsum(
        'ijk,ik->ij',
        transform[:, :2, :],
        delta
    )

    weights = np.c_[

        bary,
        1 - bary.sum(axis=1)

    ]

    rows = np.repeat(
        valid_indices,
        3
    )

    cols = vertices.flatten()

    vals = weights.flatten()

    W_linear = csr_matrix(

        (vals, (rows, cols)),
        shape=(
            target_points.shape[0],
            source_points.shape[0]
        )

    )

    return W_linear


# ============================================================
# ----------- BUILD NEAREST INTERPOLATION MATRIX -------------
# ============================================================

def build_nearest_interpolation_matrix(
    source_points,
    target_x,
    target_y
):
    
    """
   Builds the sparse nearest-neighbor interpolation operator
   used as a fallback outside the Delaunay triangulation domain.

   This operator is only applied in regions where the
   linear interpolation produces NaN values.
   """

    target_points = np.column_stack((

        target_x.flatten(order='F'),
        target_y.flatten(order='F')

    ))

    tree = cKDTree(source_points)

    _, nearest_ids = tree.query(target_points)

    rows = np.arange(target_points.shape[0])

    cols = nearest_ids

    vals = np.ones(target_points.shape[0])

    W_nearest = csr_matrix(

        (vals, (rows, cols)),
        shape=(
            target_points.shape[0],
            source_points.shape[0]
        )

    )

    return W_nearest


# ============================================================
# -------------------- APPLY INTERPOLATION -------------------
# ============================================================

def apply_interpolation(
    W_linear,
    W_nearest,
    field_unordered,
    ni,
    nj
):
    
    """
   Interpolates an unordered flow field onto the structured mesh,

       q_ordered = W q_unordered

   using the sparse linear interpolation operator.

   A nearest-neighbor correction is applied in regions
   where the linear interpolation is undefined.
   """

    field_flat = (
        W_linear @ field_unordered.flatten(order='F')
    )

    field_nearest = (
        W_nearest @ field_unordered.flatten(order='F')
    )

    mask = np.isnan(field_flat)

    if np.any(mask):

        field_flat[mask] = field_nearest[mask]

    return field_flat.reshape(
        (ni, nj),
        order='F'
    )