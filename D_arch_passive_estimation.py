def D_arch_passive_estimation(UYS_passive, rhob_passive, theta):

    # UYS was estimated the passive state at the outlet of the conical part of the hopper.
    # This UYS corresponds to the MPS, which is estimated using radial stress theory
    g = 9.8
    H = 2 + theta / 60
    D_arching_passive = H * UYS_passive[-1] / (rhob_passive[-1] * g)                                                    # arching diameter (m) in the active mode: calculation based on Eq. (2) and Eq. (3) of Leung et al paper
    D_arching_passive = D_arching_passive * 1000                                                                        # arching diameter (mm) in the active mode

    return D_arching_passive