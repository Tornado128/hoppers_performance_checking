def D_arch_active_estimation(UYS,rhob,theta):

    # UYS was estimated the active state at the outlet of the conical part of the hopper.
    # The corresponding MPS was estimated using Motzkus equations
    g = 9.8
    H = 2 + theta/60
    D_arching_active = H*UYS[-1]/(rhob[-1]*g)                                                                           # arching diameter (m) in the active mode: calculation based on Eq. (2) and Eq. (3) of Leung et al paper
    D_arching_active = D_arching_active * 1000                                                                          # arching diameter (mm) in the active mode

    return D_arching_active