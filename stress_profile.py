## This function estimates the stress profile in the hopper in active (initial discharging and filling) mode
def stress_profile(HEIGHT, RADIUS, X, Y):

    import numpy as np
    from Janssen_Equation import Janssen_Equation
    from Motzkus_Equation import Motzkus_Equation

    XX = [x for x, y in zip(X, Y) if y < HEIGHT]                            # We eliminate the x points above the powder level
    YY = [y for y in Y if y < HEIGHT]                                       # We eliminate the y points above the powder level

    YY = [HEIGHT] + YY                                                      # Height of the powder in the hopper
    XX = [RADIUS] + XX                                                      # x-location corresponds to the height of the powder in the hopper

    N = 10000                                                               # number of mesh in the vertical section
    number = len(XX) - 1                                                    # number of sections in hopper (a hopper may have several cylinderical and cone parts)
    sigmav = np.zeros(number*N)                                             # We initialized sigma0 (vertical stress) for the top of the powder as 0.
    sigma1 = np.zeros(number*N)                                             # We initialized sigma0 (major principal stress) for the top of the powder as 0.
    z = np.zeros(number*N)                                                  # (Hieght) From the top to the bottom direction (m)
    numerator = 0                                                           # numerator, which is equal to the total number of elements in z direction

    # This loop goes through all the cylinderical and cone parts of the hopper.
    # If the section is cylinderical, it uses Janssen equation to obtain vertical stress distribution.
    # if the section is cone shaped, it uses Motzkus equarion to obtain vertical stress distribution.
    for i in range(number):
        X1 = XX[i]
        X2 = XX[i + 1]
        Y1 = YY[i]
        Y2 = YY[i + 1]
        for j in range(N):
            z[numerator] = Y1 - (Y1 - Y2) * j / N
            numerator = numerator + 1
        if (X1 / X2 == 1):                                                  # If the ratio of the points remains as 1, it means that we are dealing with a cylinderical part and we will use Janssen equation to estimate the consolidation stress vs height
            sigmav_init = sigmav[i*N-i]
            [sigmav_out, sigma1_out]=Janssen_Equation(X1,X2,Y1,Y2,N,sigmav_init)
            print(sigmav_out[-1])
            sigmav[i*N:(i+1)*N] = sigmav_out[:N]
            sigma1[i*N:(i+1)*N] = sigma1_out[:N]
        else:                                                               # We are in the cone part of the hopper
            sigmav_init = sigmav[i*N-i]
            [sigmav_out, sigma1_out]=Motzkus_Equation(X1,X2,Y1,Y2,N,sigmav_init)
            print(sigmav_out[-1])
            sigmav[i*N:(i+1)*N] = sigmav_out[:N]
            sigma1[i*N:(i+1)*N] = sigma1_out[:N]
    return z, sigmav, sigma1


