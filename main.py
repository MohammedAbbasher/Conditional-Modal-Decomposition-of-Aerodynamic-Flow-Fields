# -*- coding: utf-8 -*-

import glob
import numpy as np

from load_mean_flow import load_mean_flow
from calculate_phi_cl import calculate_phi_cl
from calculate_acl_and_phicd import calculate_acl_and_phicd
from calculate_phi_cd import calculate_phi_cd
from calculate_acd import calculate_acd


def main():

    # Structured mesh dimensions
    # ni and nj correspond to the maximum grid indices
    # in the i and j directions of the Pointwise mesh   
    ni = 601
    nj = 301

    ni_ansys = ni - 1
    nj_ansys = nj - 1

    # Where the periodic cycle starts
    start_index = 24004

    # Sampling step between snapshots
    # Due to memory limitations, the flow variables were saved
    # every 5 timesteps relative to the CL and CD signals
    sampling_step = 5

    (

        W_linear,
        W_nearest,

        p_unordered,
        u_unordered,
        v_unordered,

        p_mean,
        u_mean,
        v_mean

    ) = load_mean_flow(
        ni,
        nj,
        ni_ansys,
        nj_ansys
    )

    cl1 = np.loadtxt(

        'cl18.out',
        delimiter=',',
        skiprows=1

    )

    cl = cl1[start_index::sampling_step, 1]

    cl_mean = np.mean(cl)

    print('Cl mean = ', cl_mean)

    cd1 = np.loadtxt(

        'cd18.out',
        delimiter=',',
        skiprows=1

    )

    cd = cd1[start_index::sampling_step, 1]

    cd_mean = np.mean(cd)

    print('Cd mean = ', cd_mean)

    files = sorted(
        glob.glob('p18-*')
    )

    print('Number of snapshots = ', len(files))

    (

        phi_cl,

        p_cl,
        u_cl,
        v_cl

    ) = calculate_phi_cl(

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

    )

    (

        acl,

        p_cd1,
        u_cd1,
        v_cd1,

        n_cd

    ) = calculate_acl_and_phicd(

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

    )

    (

        phi_cd,

        p_cd,
        u_cd,
        v_cd

    ) = calculate_phi_cd(

        p_cd1,
        u_cd1,
        v_cd1,

        n_cd,

        ni,
        nj

    )

    acd = calculate_acd(

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

    )

    print('\n================ Diagnostics ================')

    print('acl min = ', np.nanmin(acl))
    print('acl max = ', np.nanmax(acl))

    print('acd min = ', np.nanmin(acd))
    print('acd max = ', np.nanmax(acd))

    print('NaNs in acl = ', np.isnan(acl).sum())
    print('NaNs in acd = ', np.isnan(acd).sum())


if __name__ == "__main__":

    main()