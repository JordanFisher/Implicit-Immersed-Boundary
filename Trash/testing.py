import sys, os
root = sys.path[0]
sys.path.append(os.path.join(root, "IB_c\Release"))
import IB_c

import numpy as numpy

from Lookup import *
from FluidSolve import *
from Delta import *
from Fiber_2D import *
import pylab as m
from MatrixUtility import *
from copy import deepcopy
from time import *

from CG import CG

DeltaType = 0

def ChildN(Panel):
    for i in range(8):
        print "Child(" + str(i) + "), N =", Panel.Child(i).N

# Setup for saving the simulation output
SaveDataFlag, SaveDir, SavePeriod = True, "Sphere_Implicit", 1 # (T / dt) / 120
if SaveDataFlag:
    import SaveData as SaveData
    Save = SaveData.SaveDataControl(root, SaveDir)

# Simulation parameters
N, Nb, s, dt, T = 128, int(1**.5*64*2.**.5), 10e11, .002, .25
##N, Nb, s, dt, T = 32, 32, 10e11, .002, .25
h = 1. / N

# Plate parameters
width = .5
bottom = left = .25
y = .5

# Initialize the fiber
##X = Fiber_2D((Nb+1)**2,s)
##X.MakeSquare(Nb+1, [left,bottom,y],[left+width,bottom,y],[left+width,bottom+width,y], lambda x,y: 0*sin(x*pi)*sin(y*pi))

##y = .25
##X = Fiber_2D((Nb+1)**2,s)
##X.MakeSquare(Nb+1, [left,bottom,y],[left+width,bottom,y],[left+width,bottom+width,y], lambda x,y: .5*sin(x*pi)*sin(y*pi))

Points, Triangles = MakeSphere2(int(6000**.5),[.5,.5,.5], .2)
X = Fiber_2D(Points = Points)
X.Tether

##X = Fiber_2D((Nb+1+1) * (Nb+1)/2,s)
##MakeTriangle(Nb, [.2,.1,.5],[.8,.4,.5],[.5,.6,.5], X.X)
##Nb = 20

##Points, Triangles = MakeSphere([.5,.5,.5], .2)
##X = Fiber_2D(Points = Points)

##X = Fiber_2D(2, s)
##e = 0
##X.X[0] = [8 + e, 8 + e, 16 + e]
####X.X[0] = [15, 15, 15]
##e = 0
##X.X[1] = [12 + e, 8 + e, 16 + e]
####X.X[1] = [15.75, 15.75, 15.75]
##X.X /= float(N)

##X.X = floor(X.X * N) / N


# Initialize the fluid solver
u = Fluid(dt = dt, N = [N,N,N], dims = [1., 1., 1.], DeltaType = DeltaType)

# Induced flow parameters
flow_source = ones((N,N,N,3), float64)
Current = array([0, 0, 1], float64)

# Tree code parameters
decomp_root = os.path.join(root, "Decomposition")
Terms = 10
Min_U_Width = 2
MinPointsToHaveChildren = 10

# Initialize the fluid Green function lookup table
UField = GetLookupTable(u)

# Initialize the lookup tables needed by the tree code
import DecompositionGroup as G
import Panel as _Panel
DecompGroup = G.DecompositionGroup(decomp_root, u, Terms, Min_U_Width)


def FastFluidSolve():
    t = time()
    TreeRoot.CalcSeriesTerms(X.F)
    print "Series:", time() - t
    TreeRoot.FastEval(holdX, X.F, X.U)
    X.U *= X.hb2

def CG_FastEval(x):
    X.X = x
    X.F = -s * X.X
    FastFluidSolve()
    return x - dt * X.U

def FastEval(x):
    X.X = x
    X.CalcForceDistribution()
    FastFluidSolve()
    return x - dt * X.U

def FluidSolve():
    u.u *= 0
    FiberToGrid (u.N, u.h, X.Nb, X.hb2, holdX, X.F, u.f, DeltaType)
    u.FluidSolve(u.f)
    GridToFiber (u.N, u.h, X.Nb, X.hb2, holdX, u.u, X.U, DeltaType)

def CG_Eval(x):
    X.X = x
    X.F = -s * X.X    
    FluidSolve()
    return x - dt * X.U

def Eval(x):
    X.X = x
    X.CalcForceDistribution()
    FluidSolve()
    return x - dt * X.U


##X.F[0,:] = [1,0,0]
##for i in range(X.Nb):
##    X.F[i] = [cos(i), sin(i), cos(sin(i))]
X.F[:] = 1

holdX = X.X.copy()

t = time()
TreeRoot = _Panel.Panel(MinPointsToHaveChildren, DecompGroup.Depth, holdX, array([0,0,0], float64), array([1,1,1], float64), 0)
print "Init =", time() - t
TreeRoot.CalcUV(Terms, DecompGroup)
print "CalcUV =", time() - t
TreeRoot.PreEval(UField, DecompGroup, holdX)
print "PreEval =", time() - t
##raise
def Time(f):
    t = time()
    for i in range(1):
        f()
    print time() - t
Time(FastFluidSolve)
z = X.U.copy()
Time(FluidSolve)
w = X.U.copy()

print z[1]
print w[1]
print abs(z-w).max()

def Interp(f, i, shift):
    shiftx = shift
    shifty = X.N * shift
    return (f[i+shiftx] + f[i-shiftx] + f[i+shifty] + f[i-shifty]) / 4
    

##raw_input("")

##TreeRoot.CalcSeriesTerms(X.F)
##print TreeRoot.EvalPoint(X.X[1], UField, DecompGroup, X.X, X.F)*X.hb
##raw_input("")

##############
##############Time_FluidSolve = 0
##############Time_Initializing = 0
##############Time_Solving = 0
##############Time_Total = 0
##############Num_CGIterations = 0
##############
##############CurStep = 0
##############t = 0
##############while t < T:
##############    print "########################################################################"
##############    print "CurStep =", CurStep, ", t =", t
##############    print "argmax(X.z) =", X.X[:,2].argmax(), ", X[argmax] =", X.X[X.X[:,2].argmax()]
##############    i = (X.N-1)/2 + (X.N-1)/2 * X.N
##############    print "X[", i, "] =", X.X[i]
##############    print "u[N/2,N/2,N/4] =", u.u[N/2,N/2,N/4]
##############    print "u[N/2,N/2,3N/4] =", u.u[N/2,N/2,3*N/4]
##############    print "u.z.max() =", u.u[:,:,:,2].max(), ", u.z.min() =", u.u[:,:,:,2].min()
##############    print
##############    print
##############    angle = 0 + pi / 4 * cos(t / T * 6 * pi)
##############    Current = 100 * array([0, sin(angle), cos(angle)])    
##############    print "Current angle:", 180 * angle / pi, "   Current:", Current
##############    print
##############
##############    # Save current fiber configuration and fluid velocity
##############    clock_t1 = time()
##############    holdX = X.X.copy()
##############    holdu = u.u.copy()
##############
##############    # Initialize the tree structure
##############    clock_t2 = time()
##############    TreeRoot = _Panel.Panel(MinPointsToHaveChildren, DecompGroup.Depth, holdX, array([0,0,0], float64), array([1,1,1], float64), 0)
##############    
##############    TreeRoot.CalcUV(Terms, DecompGroup)
##############    TreeRoot.PreEval(UField, DecompGroup, holdX)
##############    Time_Initializing += time() - clock_t2
##############    
##############    # Calculate the explicit contributions
##############    clock_t2 = time()
##############    X.CalcForceDistribution()
##############    FiberToGrid (u.N, u.h, X.Nb, X.hb2, holdX, X.F, u.f, DeltaType)
##############    u.f += Current * flow_source
##############    u.FluidSolve(u.f)
##############    GridToFiber (u.N, u.h, X.Nb, X.hb2, holdX, u.u, X.U, DeltaType)
##############
##############    b = dt * X.U
##############    
##############    Time_FluidSolve += time() - clock_t2
##############    
##############    # Solve the implicit system (I - M)X = b
##############    clock_t2 = time()
##############    dx = 0 * X.X
##############
##############    dx, reps = CG(dx, b, CG_FastEval, .0001, ShowProgress = True)
##############    Num_CGIterations += reps
##############    Time_Solving += time() - clock_t2
##############
##############    # Update the fiber position
##############    X.X = holdX + dx
##############
##############    # Evolve the fluid
##############    clock_t2 = time()
##############    u.u = holdu.copy()
##############    X.CalcForceDistribution()
##############    FiberToGrid (u.N, u.h, X.Nb, X.hb2, holdX, X.F, u.f, DeltaType)
##############    u.f += Current * flow_source
##############    u.FluidSolve(u.f)
##############    Time_FluidSolve += time() - clock_t2
##############
##############    CurStep += 1
##############    t += dt
##############
##############    Time_Total += time() - clock_t1
##############
##############    # Pring benchmarks
##############    print
##############    print "Average fluid solve cost per timestep    :", Time_FluidSolve / CurStep
##############    print "Average initialization cost per timestep :", Time_Initializing / CurStep
##############    print "Average solving cost per timestep        :", Time_Solving / CurStep
##############    print "Average total cost per timestep          :", Time_Total / CurStep
##############    print
##############    print "                      (", (Time_Total / CurStep) * (T / dt) / 3600., ") hours"
##############    print
##############    print "CG iterations:", reps, " Average:", float(Num_CGIterations) / CurStep
##############    print
##############    print
##############    print
##############
##############    if SaveDataFlag and CurStep % SavePeriod == 0:
##############        #m.imshow(u.u[N/2,:,:,2], vmin = -2, vmax = 10, interpolation='bilinear')
##############        m.imshow((u.u[N/2,:,:,1]**2 + u.u[N/2,:,:,2]**2)**.5, vmin = -2, vmax = 10, interpolation='bilinear')
##############        z = N * X.X[X.N/2::X.N,:]
##############        m.scatter(z[:,2], z[:,1], c = 'k')
##############        Save.SaveStep(u, X, CurStep, SaveFigure = True)
##############
##############print
##############print "Total simulation time:", Total_Time
