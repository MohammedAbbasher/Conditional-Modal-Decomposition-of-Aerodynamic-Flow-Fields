# ============================================================
# ----------- CALCULATE ACL + BUILD PHI_CD -------------------
# ============================================================

import numpy as np

from interpolation import apply_interpolation


def calculate_acl_and_phicd(

    files,

    phi_cl,

    cd,
    cd_mean,

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
   Computes the lift modal coefficient

       a_cl = <Phi_cl, Phi'>

   where Phi' denotes the fluctuating flow field.

   The contribution of phi_cl is removed from the flow field
   to construct the residual fluctuations used for building phi_cd.
   """

    acl = np.zeros(len(files))

    p_cd = np.zeros((ni, nj))
    u_cd = np.zeros((ni, nj))
    v_cd = np.zeros((ni, nj))

    n_cd = 0

    for count, file in enumerate(files, start=1):

        print('Merged snapshot = ', count)

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

        phi_inst = np.vstack((

            p_fluctuated,
            u_fluctuated,
            v_fluctuated

        ))

        acl[count - 1] = np.sum(
            phi_cl * phi_inst
        )

        p_residual = (

            p_fluctuated
            - phi_cl[0:ni, :] * acl[count - 1]

        )

        u_residual = (

            u_fluctuated
            - phi_cl[ni:2*ni, :] * acl[count - 1]

        )

        v_residual = (

            v_fluctuated
            - phi_cl[2*ni:3*ni, :] * acl[count - 1]

        )

        if cd[count - 1] > cd_mean:

            p_cd += p_residual
            u_cd += u_residual
            v_cd += v_residual

            n_cd += 1

    return (

        acl,

        p_cd,
        u_cd,
        v_cd,

        n_cd

    )