# ============================================================
# -------------------- CALCULATE PHI_CD ----------------------
# ============================================================

from normalize_mode import normalize_mode


def calculate_phi_cd(

    p_cd,
    u_cd,
    v_cd,
    n_cd,
    ni,
    nj

):
    
    """
  Computes the normalized drag mode phi_cd from the
  conditional average of the residual fluctuating fields.

  The residual fields correspond to the flow structures
  remaining after removing the phi_cl contribution.
  """

    p_cd /= n_cd
    u_cd /= n_cd
    v_cd /= n_cd

    phi_cd = normalize_mode(

        p_cd,
        u_cd,
        v_cd,
        ni,
        nj

    )

    return (

        phi_cd,

        p_cd,
        u_cd,
        v_cd

    )