"""Test the convergence of the treecode."""

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

N, Nb, s, dt, T = 64, int(32*2.**.5), 10e11, .002, .25
u = Fluid(dt = dt, N = [N,N,N], dims = [1., 1., 1.], DeltaType = DeltaType)

import DecompositionGroup as G
decomp_root = os.path.join(root, "Decomposition")
Terms = 10
MinPanelWidth = 2
DecompGroup = G.DecompositionGroup(decomp_root, u, Terms, MinPanelWidth)

##U = DecompGroup.UField(2)
U = DecompGroup.UField(4)
do least squares over omega_in not omega_out

def GetVec(Index):
    return array([U.AB_xx(i)[Index[0],Index[1],Index[2]] for i in range(Terms)])

def GetM(Indices):
    M = zeros((Terms,len(Indices)),float64)
    for i in range(len(Indices)):
        M[:,i] = GetVec(Indices[i])
    return M

def ShiftIndices(Indices, Shift):
    return [array(Index) + array(Shift) for Index in Indices]

import numpy.linalg as linalg
Shift = [5,5,5]
##Indices = [[0,0,0],[0,0,2],[0,2,0],[2,0,0],[2,2,2],[2,2,0],[0,2,2],[2,0,2],[1,1,1],[3,3,3],[1,2,3],[3,1,2],[2,3,1]]
Indices = [[i,j,k] for i in range(18-Shift[0]) for j in range(18-Shift[1]) for k in range(18-Shift[2])]

##M = GetM(Indices[:Terms])
M = GetM(Indices)



Indices2 = ShiftIndices(Indices, Shift)
##M2 = GetM(Indices2[:Terms])
M2 = GetM(Indices2)

A = zeros((Terms,Terms),float64)
##for i in range(Terms):
##    v = zeros(Terms,float64)
##    v[i] = 1
##
##    f = linalg.solve(M,v)
##
##    A[:,i] = dot(M2, f)

##A = dot(M2, linalg.inv(M))
A = dot(dot(M2,M.T), linalg.inv(dot(M,M.T)))

Index = [0,0,0]
for Index in Indices:
    v = dot(A, GetVec(Index))
    w = GetVec(array(Index) + Shift)
    print abs(v-w).max()

z = GetVec([35,35,37])
print dot(z,v), dot(z,w)
