"""Test the performance and accuracy of the treecode for increasing values of p."""

import sys, os
root = sys.path[0]
sys.path.append(os.path.join(root, "IB_c\Release"))
import IB_c

from numpy.linalg import eigvals
from numpy.linalg import solve

import numpy as numpy

from Lookup import *
from FluidSolve import *
from Delta import *
from Fiber_2D import *
##from Plot import *
import pylab as m
##m.show()
from MatrixUtility import *
from copy import deepcopy
from time import *

from CG import CG

DeltaType = 0

# Simulation parameters
##N, Nb, dt, s = 128, int(64)+2, .002, 10e7
N, Nb, dt, s = 128, int(64), .002, 10e7
h = 1. / N

# Plate parameters
width = .5
bottom = left = .25
y = .5

# Initialize the fiber
X = Fiber_2D((Nb+1)**2,s)
X.MakeSquare(Nb+1, [left,bottom,y],[left+width,bottom,y],[left+width,bottom+width,y], lambda x,y: 0.9*h*sin(4*x*pi)*sin(4*y*pi))

holdX = X.Tether.copy()
##holdX = floor(X.X*N)/N

# Initialize the fluid solver
u = Fluid(dt = dt, N = [N,N,N], dims = [1., 1., 1.], DeltaType = DeltaType)

# Tree code parameters
Thin = True            # Whether we are doing a thin, 2D slice of the Green's function, or the entire 3D domain
Thin_Width = 1
Thin_z = .5
if Thin: 
    decomp_root = os.path.join(root, "Thin_Decomposition")
else:
    decomp_root = os.path.join(root, "Decomposition")

Terms = 0
MinPanelWidth = 2
MinPointsToHaveChildren = 10

# Initialize the fluid Green function lookup table
UField = GetLookupTable(u, Thin_Width)


def FastFluidSolve():
    TreeRoot.CalcSeriesTerms(X.F)
    TreeRoot.FastEval(holdX, X.F, X.U)
    X.U *= X.hb2

def FastEval():
    X.CalcForceDistribution()
    FastFluidSolve()
    return dt * X.U

def FluidSolve():
    u.u *= 0
    FiberToGrid (u.N, u.h, X.Nb, X.hb2, holdX, X.F, u.f, DeltaType)
    u.FluidSolve(u.f)
    GridToFiber (u.N, u.h, X.Nb, X.hb2, holdX, u.u, X.U, DeltaType)

def Eval():
    X.CalcForceDistribution()
    FluidSolve()
    return dt * X.U

import DecompositionGroup as G
import Panel as _Panel

FluidSolve_dX = Eval()

TotalTerms = 50
Iterations = 10
InitializationTime, PreEvalTime, EvalTime = [], [], []
MaxDifference, L2Difference, AverageRelativeError = [], [], []

# Initialize the lookup tables needed by the tree code
DecompGroup = G.DecompositionGroup(decomp_root, N, TotalTerms, MinPanelWidth, Thin, Thin_Width, Thin_z)

while Terms < TotalTerms:
    print "Terms =", Terms    
    
    # Time the preevaluation and initialization of the tree structure
    clock_t = time()
    for i in range(Iterations/2):    
        TreeRoot = _Panel.Panel(MinPointsToHaveChildren, DecompGroup.Depth, holdX, array([0,0,0], float64), array([1,1,1], float64), 0)
        TreeRoot.CalcUV(Terms, DecompGroup)
        TreeRoot.PreEval(UField, DecompGroup, holdX)
    PreEvalTime.append((time() - clock_t) / (Iterations/2))

    # Time the tree code evaluation    
    clock_t = time()
    for i in range(Iterations):
        TreeRoot.CalcSeriesTerms(X.F)
        TreeRoot.FastEval(holdX, X.F, X.U)
    EvalTime.append((time() - clock_t) / Iterations)

    # Check the accuracy of the tree code
    TreeEval_dX = FastEval()
    dif = TreeEval_dX - FluidSolve_dX
    MaxDifference.append(abs(dif).max())
    L2Difference.append((dif**2).sum()**.5/X.Nb)
    AverageRelativeError.append(abs(dif/FluidSolve_dX).sum() / X.Nb)
                        

    Terms += 1

print "Initialization time"
print InitializationTime
print
print "Preevaluation time"
print PreEvalTime
print
print "Evaluation time"
print EvalTime
print
print
print
print "Max difference"
print MaxDifference
print
print "L2Difference"
print L2Difference
