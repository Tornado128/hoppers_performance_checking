## This function estimates the stress profile in the hopper in active (initial discharging and filling) mode
#def stress_in_active_state(a, b, c, X, Y, UPPER, LOWER, HEIGHT):

a = [1.46820423e+01,  7.86672933e+03,  7.93957291e-01,  2.00130185e-03, -2.52808810e+01]
b = [3.27446958e-01, 1.38287317e-04, 3.76692531e+00, 2.90532160e+01, 8.84732787e-02]
c = [125.25240252, -7812.46747149, 0, 0, 90]

Y = [0.525, 0.325, 0]
X = [0.2845, 0.2845, 0.2]

HEIGHT = 0.43720000000000003
RADIUS = 0.2845

XX = [x for x, y in zip(X, Y) if y <= HEIGHT]
YY = [y for y in Y if y <= HEIGHT]
#print(XX)
#print(YY)

YY = [HEIGHT] + YY
XX = [RADIUS] + XX
print(YY)
print(XX)

# 0 -> rhob
# 1 -> PHIE
# 2 -> FC
# 3 -> PHILIN
# 4 -> wall friction

sigma0 = 0                                                              # We initialized sigma0 (consolidation stress) for the top of the powder as 0.
number = len(XX) - 1                                                    # number of sections in hopper (a hopper may have several cylinderical and cone parts)
for i in range(number):
    X1 = XX[i]
    X2 = XX[i + 1]
    Y1 = YY[i]
    Y2 = YY[i + 1]
    if (X1 / X2 == 1):                                        # If the ratio of the points remains as 1, it means that we are dealing with a cylinderical part and we will use Janssen equation to estimate the consolidation stress vs height
        Janssen_Equation(sigma0,X1,X2,Y1,Y2,a,b,c)


