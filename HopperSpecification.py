# This constructor defines several axisymmetric hoppers with their name, and their 2D sketch in the form of x and y.
# All these hoppers are conical
class ConicalHoppers:
    def __init__(self, name, z, x):
        self.name = name
        self.z = z
        self.x = x

n = 7               # number of hoppers that we have information about them
r =[0]*n            # defining a list with n elements of zero
r[0] = ConicalHoppers("1000L_IBC",[1.4, 0.883, 0], [0.6098, 0.6098, 0.1])
r[1] = ConicalHoppers("100L_IBC",[0.525, 0.325, 0], [0.2845, 0.2845, 0.2])
r[2] = ConicalHoppers("Piccola",[0.385, 0.195, 0.05, 0],[0.059, 0.059, 0.015, 0.015])
r[3] = ConicalHoppers("Courtoy",[0.49, 0.28, 0.028, 0], [0.0955, 0.0955, 0.043, 0.043])
r[4] = ConicalHoppers("Kosrsch_XL100_poor_flow",[0.356, 0.326, 0.03, 0], [0.099, 0.099, 0.037, 0.037])
r[5] = ConicalHoppers("Korsch_XL100_gravity",[0.518, 0.478, 0.313, 0.125, 0],[0.125, 0.125, 0.025, 0.025, 0.018])
r[6] = ConicalHoppers("300L_IBC", [0.923, 0.643, 0.08, 0],[0.425, 0.425, 0.1, 0.1])


