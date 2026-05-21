# ============================================================
# -------------------- CALCULATE PHI_CL ----------------------
# ============================================================

import numpy as np

from normalize_mode import normalize_mode
from interpolation import apply_interpolation


def calculate_phi_cl(

    files,

    cl,
    cl_mean,

    W_linear,
    W_nearest,

    p_unordered,
    u_unordered,
    v_unordered,

    ni,
    nj,
    ni_ansys,
    nj_ansys

):
    
    
    """
    Computes the conditional lift mode phi_cl from snapshots satisfying

        CL > mean(CL)

    using conditional averaging of the fluctuating flow fields.

    The resulting mode represents the dominant flow structures
    associated with high-lift events.
    """

    p_cl = np.zeros((ni, nj))
    u_cl = np.zeros((ni, nj))
    v_cl = np.zeros((ni, nj))

    n_cl = 0

    for count, file in enumerate(files, start=1):

        print('Conditional CL snapshot = ', count)

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

        if cl[count - 1] > cl_mean:

            p_cl[0:ni_ansys, 0:nj_ansys] += (
                p - p_unordered
            )

            u_cl[0:ni_ansys, 0:nj_ansys] += (
                u - u_unordered
            )

            v_cl[0:ni_ansys, 0:nj_ansys] += (
                v - v_unordered
            )

            n_cl += 1

    print('Number of CL samples = ', n_cl)

    p_cl /= n_cl
    u_cl /= n_cl
    v_cl /= n_cl

    p_cl = apply_interpolation(
        W_linear,
        W_nearest,
        p_cl[0:ni_ansys, 0:nj_ansys],
        ni,
        nj
    )

    u_cl = apply_interpolation(
        W_linear,
        W_nearest,
        u_cl[0:ni_ansys, 0:nj_ansys],
        ni,
        nj
    )

    v_cl = apply_interpolation(
        W_linear,
        W_nearest,
        v_cl[0:ni_ansys, 0:nj_ansys],
        ni,
        nj
    )

    phi_cl = normalize_mode(
        p_cl,
        u_cl,
        v_cl,
        ni,
        nj
    )

    return (

        phi_cl,

        p_cl,
        u_cl,
        v_cl

    )