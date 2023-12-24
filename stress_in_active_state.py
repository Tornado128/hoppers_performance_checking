## This function estimates the stress profile in the hopper in active (initial discharging and filling) mode
#def stress_in_active_state(a, b, c, X, Y, UPPER, LOWER, HEIGHT):

a = [1.46820423e+01, 7.86672933e+03, 7.93957291e-01, 2.00130185e-03, 8.27337635e-01]
b = [3.27446958e-01, 1.38287317e-04, 3.76692531e+00, 2.90532160e+01, 9.51315496e+01]
c = [125.25240252, -7812.46747149, 0, 0, 0]

Y = [0.525, 0.325, 0]
X = [0.2845, 0.2845, 0.2]

HEIGHT = 0.31
UPPER = 0.2845
LOWER = 0.2845

XX = [x for x, y in zip(X, Y) if y <= HEIGHT]
YY = [y for y in Y if y <= HEIGHT]
#print(XX)
#print(YY)

YY = [HEIGHT] + YY
print(YY)
# 0 -> rhob
# 1 -> PHIE
# 2 -> FC
# 3 -> PHILIN
# 4 -> TAU (wall friction)
