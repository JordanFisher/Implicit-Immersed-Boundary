"""Test the accuracy of the matrix approximation against a direct fluid solve.
"""

from numpy import zeros, float64, floor, array, dot

import IB

from Fluid import Fluid
from Ellipse import Ellipse
from MatrixUtility import Stack, Unstack

# Simulation parameters
N = 128
RestrictToEulerian = False # If True we restrict the fiber points to Eulerian intersections. When done the matrix should give EXACT results.

# Create the fiber in the shape of an ellipse
X = Ellipse(2 * N, s = 10e6, c = [.5,.5], p1 = [.8,.5], p2 = [.5,.6])
if RestrictToEulerian:
    X.X = floor(X.X * N) / N 

# Initialize the fluid solver
u = Fluid(dt = 1. / N, N = [N,N], dims = [1., 1.])


# Calculate dX/dt (stored in X.U) using a fluid solve
X.CalcForceDistribution()
X.ToGrid(u)
u.UpdateFluid(u.f)
X.FromGrid(u)


# Calculate dX/dt (stored in U) using the fluid matrix
# First construt the fluid matrix
M = X.ConstructFluidMatrix(u)

# Then calculate U = M * F, where F is the fiber force
U = Unstack(dot(M, Stack(X.F)), 2)

# Calculate difference between X.U and U (the induced flow calculated via a fluid solve and via the matrix)
print abs(X.U - U).max() / abs(X.U).max()
