"""Test the performance and accuracy of the treecode for increasing values of Nb."""


Broken code! Not up-to-date!
Can we merge this with the other test?


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
N, Max_Nb, dt, s = 128, 100, .002, 10e7
h = 1. / N

# Plate parameters
width = .5
bottom = left = .25
y = .5


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

Terms = 10
MinPanelWidth = 2
MinPointsToHaveChildren = 10

# Initialize the fluid Green function lookup table
UField = GetLookupTable(u, Thin_Width)


def FastFluidSolve():
    TreeRoot.CalcSeriesTerms(X.F)
    TreeRoot.FastEval(X.X, X.F, X.U)
    X.U *= X.hb2

def FastEval():
    X.CalcForceDistribution()
    FastFluidSolve()
    return dt * X.U

import DecompositionGroup as G
import Panel as _Panel

Iterations = 10
InitializationTime, PreEvalTime, EvalTime = [], [], []

# Initialize the lookup tables needed by the tree code
DecompGroup = G.DecompositionGroup(decomp_root, N, Terms, MinPanelWidth, Thin, Thin_Width, Thin_z)

Nb = 90
while Nb < Max_Nb:
    print "Nb =", Nb

    # Initialize the fiber
    X = Fiber_2D(Nb**2,s)
    X.MakeSquare(Nb, [left,bottom,y],[left+width,bottom,y],[left+width,bottom+width,y], lambda x,y: 0.9*h*sin(4*x*pi)*sin(4*y*pi))
        
##    # Time the initialization of the tree structure
##    clock_t = time()
##    for i in range(Iterations):    
##        TreeRoot = _Panel.Panel(MinPointsToHaveChildren, DecompGroup.Depth, X.X, array([0,0,0], float64), array([1,1,1], float64), 0)
##    InitializationTime.append((time() - clock_t) / Iterations)
##
##    # Time the preevaluation
##    clock_t = time()
##    for i in range(Iterations):
##        TreeRoot.CalcUV(Terms, DecompGroup)
##        TreeRoot.PreEval(UField, DecompGroup, X.X)
##    PreEvalTime.append((time() - clock_t) / Iterations)

    # Time the initialization of the tree structure and pre-eval
    clock_t = time()
    for i in range(Iterations/2):
        TreeRoot = _Panel.Panel(MinPointsToHaveChildren, DecompGroup.Depth, X.X, array([0,0,0], float64), array([1,1,1], float64), 0)
        TreeRoot.CalcUV(Terms, DecompGroup)
        TreeRoot.PreEval(UField, DecompGroup, X.X)
    PreEvalTime.append((time() - clock_t) / (Iterations/2))

    # Time the tree code evaluation    
    clock_t = time()
    for i in range(Iterations):
        TreeRoot.CalcSeriesTerms(X.F)
        TreeRoot.FastEval(X.X, X.F, X.U)
    EvalTime.append((time() - clock_t) / Iterations)

    Nb += 1

##print InitializationTime
##print
print PreEvalTime
print
print EvalTime
