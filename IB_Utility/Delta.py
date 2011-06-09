"""Implements the Spreading and Interpolation operators of the Immersed Boundary Method."""

import sys
import os as os
sys.path.append(os.path.abspath("IB_c\Release"))
import IB_c
from numpy import *

def GetDeltaFunction(h, DeltaType):
	"""Return a discrete delta function.
	
	DeltaType is an integer specifying the formula to use
	DeltaType == 0:
		Delta(r) = 
			/ (1 + cos(pi * x / 2) / (4 * h) for |x| < 2
			\ 0									otherwise
	DeltaType == 1:
		Delta(r) =
			/  .125 * (3 - 2 * |x| + sqrt(1 + 4 * |x| - 4x**2)) / h for 0 <= |x| < 1
			|  .125 * (5 - 2 * |x| - sqrt(-7 + 12 * |x| - 4x**2)) / h for 1 < |x| <= 2
			\  0									otherwise
	where x = r / h"""

	if DeltaType == 0:
		def Delta(r):
			if abs(r) < 2 * h:
				return (1. + cos(pi * r / (2*h))) / (4. * h)
			else:
				return 0
	else:
		def Delta(r):
			x = r / h
			absx = abs(x)
			if absx <= 2:
				if absx <= 1:
					return .125 * (3. - 2 * absx + sqrt(1. + 4. * absx - 4. * x * x)) / h
				else:
					return .125 * (5. - 2 * absx - sqrt(-7. + 12. * absx - 4 * x * x)) / h
			else:
				return 0
				
	return Delta
	
	
def Delta(h, r, DeltaType):
	"""Compute the 1D discrete delta approximation at h.
	
	DeltaType is an integer specifying the formula to use
	DeltaType == 0:
		(1 + cos(pi * x / 2) / (4 * h) for |x| < 2
		0									otherwise
	DeltaType == 1:
		.125 * (3 - 2 * |x| + sqrt(1 + 4 * |x| - 4x**2)) / h for 0 <= |x| < 1
		.125 * (5 - 2 * |x| - sqrt(-7 + 12 * |x| - 4x**2)) / h for 1 < |x| <= 2
		0									otherwise
	where x = r / h"""
		
	if DeltaType == 0:
		if abs(r) < 2 * h:
			return (1. + cos(pi * r / (2*h))) / (4. * h)
		else:
			return 0
	else:
		x = r / h
		absx = abs(x)
		if absx <= 2:
			if absx <= 1:
				return .125 * (3. - 2 * absx + sqrt(1. + 4. * absx - 4. * x * x)) / h
			else:
				return .125 * (5. - 2 * absx - sqrt(-7. + 12. * absx - 4 * x * x)) / h
		else:
			return 0

def FiberToGrid (N, h, Nb, hb, X, F, f, DeltaType = 0, UseNativePython = False):
	"""Spread a distribution F on X to the Eulerian grid f."""

	Dim = len(N)
	
	if Dim == 2 or Dim == 3:
		return FiberToGrid_2D(N, h, Nb, hb, X, F, f, DeltaType, UseNativePython)
	elif Dim == 3:
		return FiberToGrid_3D(N, h, Nb, hb, X, F, f, DeltaType, UseNativePython)
	else:
		raise Exception("N must be of length 2 or 3.")
	
def FiberToGrid_2D (N, h, Nb, hb, X, F, f, DeltaType = 0, UseNativePython = False):
	"""Spread a distribution F on X to the Eulerian grid f."""
	
	f *= 0

	if UseNativePython:
		# Loop through fiber points
		for i in range(Nb):
			# Modify force for nearby points
			jMin = int(X[i][0] / h[0] - 2.)
			jMax = jMin + 5
			kMin = int(X[i][1] / h[1] - 2.)
			kMax = kMin + 5
		   
			for j in range(jMin, jMax + 1):
				_delt = Delta(h[0], j * h[0] - X[i][0], DeltaType) * hb
				for k in range(kMin, kMax + 1):
					delt = _delt * Delta(h[1], k * h[1] - X[i][1], DeltaType)
					f[j%N[0],k%N[1]] += F[i] * delt
	else:
		IB_c.FiberToGrid(N, h, Nb, hb, X, F, f, DeltaType)
		
	return f
	
def FiberToGrid_3D (N, h, Nb, hb, X, F, f, DeltaType = 0, UseNativePython = False):
	"""Spread a distribution F on X to the Eulerian grid f."""
	
	f *= 0

	if UseNativePython:
		# Loop through fiber points
		for i in range(Nb):
			# Modify force for nearby points
			jMin = int(X[i][0] / h[0] - 2.)
			jMax = jMin + 5
			kMin = int(X[i][1] / h[1] - 2.)
			kMax = kMin + 5
			lMin = int(X[i][2] / h[2] - 2.)
			lMax = lMin + 5
		   
			for j in range(jMin, jMax + 1):
				__delt = Delta(h[0], j * h[0] - X[i][0], DeltaType) * hb
				for k in range(kMin, kMax + 1):
					_delt = __delt * Delta(h[1], k * h[1] - X[i][1], DeltaType)
					for l in range(lMin, lMax + 1):
						delt = _delt * Delta(h[2], l * h[2] - X[i][2], DeltaType)
						f[j%N[0],k%N[1],l%N[2]] += F[i] * delt
	else:
		IB_c.FiberToGrid(N, h, Nb, hb, X, F, f, DeltaType)
		
	return f


def GridToFiber (N, h, Nb, hb, X, u, U, DeltaType = 0, UseNativePython = False):
	"""Interpolate a vector-field u to the points in X.
	
	The result is stored in U."""

	Dim = len(N)
	
	if Dim == 2:
		return GridToFiber_2D(N, h, Nb, hb, X, u, U, DeltaType, UseNativePython)
	elif Dim == 3:
		return GridToFiber_3D(N, h, Nb, hb, X, u, U, DeltaType, UseNativePython)
	else:
		raise Exception("N must be of length 2 or 3.")	

def GridToFiber_2D(N, h, Nb, hb, X, u, U, DeltaType = 0, UseNativePython = False):
	"""Interpolate a vector-field u to the points in X.
	
	The result is stored in U."""
	
	U *= 0

	if UseNativePython:
		# Loop through fiber points
		for i in range(Nb):
			# Modify velocity for nearby points
			jMin = int(X[i][0] / h[0] - 2.)
			jMax = jMin + 5
			kMin = int(X[i][1] / h[1] - 2.)
			kMax = kMin + 5

			for j in range(jMin, jMax + 1):
				_delt = Delta(h[0], j * h[0] - X[i][0], DeltaType) * h[0] * h[1]
				for k in range(kMin, kMax + 1):
					delt = _delt * Delta(h[1], k * h[1] - X[i][1], DeltaType)
					U[i] += u[j%N[0],k%N[1]] * delt
	else:
		IB_c.GridToFiber(N, h, Nb, hb, X, u, U, DeltaType)
		
	return U

def GridToFiber_3D(N, h, Nb, hb, X, u, U, DeltaType = 0, UseNativePython = False):
	"""Interpolate a vector-field u to the points in X.
	
	The result is stored in U."""
	
	U *= 0

	if UseNativePython:
		# Loop through fiber points
		for i in range(Nb):
			# Modify velocity for nearby points
			jMin = int(X[i][0] / h[0] - 2.)
			jMax = jMin + 5
			kMin = int(X[i][1] / h[1] - 2.)
			kMax = kMin + 5
			lMin = int(X[i][2] / h[2] - 2.)
			lMax = lMin + 5

			for j in range(jMin, jMax + 1):
				__delt = Delta(h[0], j * h[0] - X[i][0], DeltaType) * h[0] * h[1] * h[2]
				for k in range(kMin, kMax + 1):
					_delt = __delt * Delta(h[1], k * h[1] - X[i][1], DeltaType)
					for l in range(lMin, lMax + 1):
						delt = _delt * Delta(h[2], l * h[2] - X[i][2], DeltaType)
						U[i] += u[j%N[0],k%N[1],l%N[2]] * delt
	else:
		IB_c.GridToFiber(N, h, Nb, hb, X, u, U, DeltaType)
		
	return U