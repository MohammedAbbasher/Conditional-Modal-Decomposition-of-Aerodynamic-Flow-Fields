# -*- coding: utf-8 -*-

import numpy as np


def normalize_mode(p, u, v, ni, nj):

    phi1 = np.vstack((p, u, v))

    phi_phi = np.sqrt(
        np.sum(
            (phi1.T @ phi1),
            axis=1
        )
    )

    phi = np.zeros_like(phi1)

    for j in range(nj):

        if phi_phi[j] != 0:

            phi[:, j] = (
                phi1[:, j] / phi_phi[j]
            )

    return phi