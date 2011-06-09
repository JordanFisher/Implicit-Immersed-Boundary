import IB_c

import numpy as numpy
from numpy import array, zeros, ones, float64, complex64, fft, real

from Lookup import GetLookupTable
from SaveData import ensure_dir, ensure_filedir

import os
from os import path

DeltaType = 0
	
PowerIterations = 10	# How many iterations are used to approximate each singular value

def GetDecompPath():
	"""Return the full path to the folder holding the Decomposition."""
	
	return path.join(path.dirname(path.dirname(__file__)), "Decomposition")
		
def GetDecomposition (fluid, Terms = 10, MinPanelWidth = 1):
	"""Get the far field decomposition associated with a fluid.
	
	If it does not exist, prompt the user if it should be created."""
	
	DecompPath = GetDecompPath()
	
	import DecompositionGroup as G
	
	# Try loading the decomposition
	try:
		DecompGroup = G.DecompositionGroup(DecompPath, fluid, Terms, MinPanelWidth)
	except Exception, (msg):
		folder = G.DecompFolderRoot(DecompPath, fluid)
		print "Decomposition group at\n	", folder, "\ncouldn't be loaded.\n"
		print msg

		response = raw_input("Would you like to create it? (y/n)")
	
		if response[0] == 'n' or response[0] == 'N':
			raise Exception("Implicit simulation requires a decomposition group")

		print
		CreateDecomposition(fluid, Terms, MinPanelWidth)
		print "Done."

		DecompGroup = G.DecompositionGroup(DecompPath, fluid, Terms, MinPanelWidth)

	return DecompGroup

def Shift(f, x, y, z):
	"""Take a periodic scalar field f and translate it by (x,y,z), where (x,y,z) is a set of integer indices."""
	
	w = f * 0

	w[x:, y:, z:] = f[:-x, :-y, :-z]
	w[x:, :y, z:] = f[:-x, -y:, :-z]
	w[:x, y:, z:] = f[-x:, :-y, :-z]
	w[:x, :y, z:] = f[-x:, -y:, :-z]
	w[x:, y:, :z] = f[:-x, :-y, -z:]
	w[x:, :y, :z] = f[:-x, -y:, -z:]
	w[:x, y:, :z] = f[-x:, :-y, -z:]
	w[:x, :y, :z] = f[-x:, -y:, -z:]

	return w

def Convolve(xhat, y):
	"""Convolve the scalar field y with x. xhat is the Fourier transform of x."""
	
	yhat = fft.fftn(y)
	return real(fft.ifftn(xhat * yhat))

			
def CreateDecomposition(u, NumTerms = 50, MinPanelWidth = 1):
	N, dt = u.N[0], u.dt

	# Get the path for where to save the decomposition
	DecompPath = GetDecompPath()

	# Get the lookup table for the Greens function of the Spread/Interpolated Stokes flow
	Greens = GetLookupTable(u)


	# The Greens function is a 3x3 matrix field, hence we must calculate 9 expansions, one for each component.
	# By symmetry only 3 of these expansions are unique.
	for __j in range(1):
		for __k in range(2):
			# Get the (__j,__k) component of the Greens function and call it Phi, and its Fourier transform PhiHat
			Phi = Greens[...,__j,__k]
			PhiHat = real(fft.fftn(Phi))

			# The largest sized panel we care about has width N / 4, because larger panels are never well separated from anything.
			PanelWidth = N / 4
			
			# We loop through all panel sizes, starting with N / 4, then N / 8, then N / 16, and so on, until a preset stopping size is passed.
			while PanelWidth > MinPanelWidth:
				print "__j ==", __j
				print "__k ==", __k
				print "PanelWidth ==", PanelWidth
				
				# Two lists of expansion coefficients, one for the incoming and one for the outgoing.
				U, V = [], []
				
				# Fill the lists with arbitrary initial guesses.
				for l in range(NumTerms):
					U.append(ones((N,N,N)))
					V.append(ones((N,N,N)))
					for j1 in range(N):
						for j2 in range(N):
							for j3 in range(N):
								V[l][j1,j2,j3] = 1

				# Define the incoming and outgoing regions.
								
				# Outgoing Region (Omega_Out)
				# [U_x1,U_x2]x[U_y1,U_y2]x[U_z1,U_z2]
				U_x1 = N / 2 - PanelWidth / 2
				U_x2 = N / 2 + PanelWidth / 2
				U_y1 = N / 2 - PanelWidth / 2
				U_y2 = N / 2 + PanelWidth / 2
				U_z1 = N / 2 - PanelWidth / 2
				U_z2 = N / 2 + PanelWidth / 2
				U_x1 -= 1
				U_x2 += 1
				U_y1 -= 1
				U_y2 += 1
				U_z1 -= 1
				U_z2 += 1

				# Incoming Region (Omega_In)
				# [V_x1,V_x2]x[V_y1,V_y2]x[V_z1,V_z2]
				V_width = 3 * PanelWidth
				V_x1 = N / 2 - V_width / 2
				V_x2 = N / 2 + V_width / 2
				V_y1 = N / 2 - V_width / 2
				V_y2 = N / 2 + V_width / 2
				V_z1 = N / 2 - V_width / 2
				V_z2 = N / 2 + V_width / 2
				V_x1 += 1
				V_x2 -= 1
				V_y1 += 1
				V_y2 -= 1
				V_z1 += 1
				V_z2 -= 1
				
				# Define the indicator functions for the incoming and outgoing regions.
				
				def Enforce_U_Restriction(U):
					"""Ensure that the tensor field U is zero outside the Outgoing Region"""
					hold = U[U_x1:U_x2+1, U_y1:U_y2+1, U_z1:U_z2+1].copy()
					U *= 0
					U[U_x1:U_x2 + 1,U_y1:U_y2 + 1,U_z1:U_z2 + 1] = hold

					return U

				def Enforce_V_Restriction(V):
					"""Ensure that the tensor field V is zero outside the Incoming Region"""
					V[V_x1:V_x2 + 1,V_y1:V_y2 + 1,V_z1:V_z2 + 1] *= 0
					V[1:N/2,:,:] *= 0
					V[:,1:N/2,:] *= 0
					V[:,:,1:N/2] *= 0

					return V

				U[0] = Enforce_U_Restriction(U[0])
				
				# A list to store the singular values
				_s = []
				
				# We wish to calculate the first few singular values, from 1 to NumTerms
				for l in range(NumTerms):
					V[l] = Enforce_V_Restriction(V[l])
					
					# We perform the Power Iteration a fixed number of times to approximate each singular value.
					for i in range(PowerIterations):
						# Calculate the product Phi' * U
						less = 0
						for _l in range(l):
							less += U[_l]*(V[_l]*V[l]).sum()
						U[l] = (Convolve(PhiHat, V[l]) - less) / (V[l]**2).sum()
						U[l] = Enforce_U_Restriction(U[l])

						# Renormalize
						norm1 = (U[l]**2).sum()**.5
						norm2 = (V[l]**2).sum()**.5
						U[l] *= (norm2 / norm1)**.5
						V[l] *= (norm1 / norm2)**.5

						# Calculate the product Phi' * V
						less = 0
						for _l in range(l):
							less += U[_l]*(V[_l]*U[l]).sum()
						V[l] = (Convolve(PhiHat, U[l]) - less) / (U[l]**2).sum()
						V[l] = Enforce_V_Restriction(V[l])

						# Renormalize
						norm1 = (U[l]**2).sum()**.5
						norm2 = (V[l]**2).sum()**.5
						U[l] *= (norm2 / norm1)**.5
						V[l] *= (norm1 / norm2)**.5
						
					print l,
						
					_s.append( ((U[l]**2).sum() * (V[l]**2).sum())**.5 )

				# Save the decomposition
				for l in range(NumTerms):
					name = path.join(DecompPath, u.Params_ToString())
					name = path.join(name, 'Panel_Width_' + str(PanelWidth))
					name = path.join(name, str(l) + '_' + str(__j) + str(__k) + '.npy')
					ensure_filedir(name)

					# We store both components U and V together in the same scalar field.
					Composite = U[l] + V[l]
					Composite = Shift(Composite, -2, -2, -2)
					numpy.save(name, Composite[N/2 - PanelWidth/2 - 2:, N/2 - PanelWidth/2 - 2:, N/2 - PanelWidth/2 - 2:])

				# Sort and print the singular values
				_s.sort()
				_s.reverse()
				print "\n\nSingular Values:\n" + str(_s) + "\n\n\n"

				# Reduce the width of Omega_Out by half
				PanelWidth /= 2


				
				
## Test code. Use these to test the accuracy of the expansions.
def Eval(j, k, Group, Layer = 2, Component = 'xx'):
	"""Evaluate the influence between two positions j and k using an expansion given by Group.
	
	j and k are 3-vectors of integer indices."""
	
	D = Group.UField(Layer)
	
	result = 0
	for i in range(D.Terms):
		if Component == 'xx':			
			result += D.AB_xx(i)[k[0],k[1],k[2]] * D.AB_xx(i)[j[0], j[1], j[2]]		
		if Component == 'xy':			
			result += D.AB_xy(i)[k[0],k[1],k[2]] * D.AB_xy(i)[j[0], j[1], j[2]]		

	return result

def EvalDirect(j, k, Greens, Component = 'xx'):
	"""Evaluate the influence between two positions j and k using a direct lookup in the Greens function lookup table.
	
	j and k are 3-vectors of integer indices."""
	
	if Component == 'xx':
		return Greens[j[0]-k[0], j[1]-k[1], j[2]-k[2], 0, 0]
	else:
		return Greens[j[0]-k[0], j[1]-k[1], j[2]-k[2], 0, 1]

### G is a DecompGroup. These two different evaluations should be close in value.
##G.Eval([0,0,0],[20,20,20],DecompGroup)
##G.EvalDirect([20,20,20],[0,0,0],Greens)