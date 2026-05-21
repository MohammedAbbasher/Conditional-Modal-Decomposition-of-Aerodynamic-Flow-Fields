# ============================================================
# ---------------------- CALCULATE ACD -----------------------
# ============================================================

import numpy as np

from interpolation import apply_interpolation


def calculate_acd(

    files,

    phi_cl,
    phi_cd,

    W_linear,
    W_nearest,

    p_mean,
    u_mean,
    v_mean,

    ni,
    nj,
    ni_ansys,
    nj_ansys

):
    
    """
  Computes the drag modal coefficient

      a_cd = <Phi_cd, Phi_res>

  where Phi_res represents the residual fluctuating field
  obtained after subtracting the phi_cl contribution.

  The coefficient quantifies the temporal evolution of the
  drag-related coherent flow structures.
  """

    acd = np.zeros(len(files))

    for count, file in enumerate(files, start=1):

        print('acd snapshot = ', count)

        struct = np.loadtxt(

            file,
            delimiter=',',
            skiprows=1

        )

        data = np.reshape(

            struct,
            (ni_ansys, nj_ansys, 6),
            order='F'

        )

        p = data[:, :, 3]
        u = data[:, :, 4]
        v = data[:, :, 5]

        p_ordered = apply_interpolation(
            W_linear,
            W_nearest,
            p,
            ni,
            nj
        )

        u_ordered = apply_interpolation(
            W_linear,
            W_nearest,
            u,
            ni,
            nj
        )

        v_ordered = apply_interpolation(
            W_linear,
            W_nearest,
            v,
            ni,
            nj
        )

        p_fluctuated = (
            p_ordered - p_mean
        )

        u_fluctuated = (
            u_ordered - u_mean
        )

        v_fluctuated = (
            v_ordered - v_mean
        )

        acl_local = np.sum(

            phi_cl * np.vstack((

                p_fluctuated,
                u_fluctuated,
                v_fluctuated

            ))

        )

        p_residual = (

            p_fluctuated
            - phi_cl[0:ni, :] * acl_local

        )

        u_residual = (

            u_fluctuated
            - phi_cl[ni:2*ni, :] * acl_local

        )

        v_residual = (

            v_fluctuated
            - phi_cl[2*ni:3*ni, :] * acl_local

        )

        phi_res = np.vstack((

            p_residual,
            u_residual,
            v_residual

        ))

        acd[count - 1] = np.sum(
            phi_cd * phi_res
        )

    return acd