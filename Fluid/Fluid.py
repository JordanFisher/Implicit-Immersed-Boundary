from numpy import mgrid, array, zeros, int32, float64, complex64, fft
import pylab as m

import Symbols
import Derivatives

from Lookup import GetLookupTable
from FieldUtility import FieldDot, FieldFFT, FieldIFFT

from IB_c import Eno2, Eno3

class ConvectionMethods:
	Centered, Upwind, ENO2, ENO3 = range(4)
	Name = { Centered : 'Centered', Upwind : 'Upwind', ENO2 : 'ENO2', ENO3 : 'ENO3' }

class PlotStyle:
	"""Different choices for plotting the fluid values."""
	
	vel, magnitude, xvel, yvel, zvel, pressure, vorticity = range(7)
	
class Fluid:
	"""Holds all the properties and data associated to a boxed fluid."""
	
	def __init__(self, dt, N, dims, rho = 1., mu = 1., DeltaType = 0, Dim = None,
	             ConvectionMethod = ConvectionMethods.Centered):
		"""dims should be a 2-vector or 3-vector holding the dimensions of the fluid box.
		N should likewise be an integer 2-vector or 3-vector which describes how fine the grid will be in each dimension."""

		if Dim == None:
			self.Dim = Dim = len(N)
		else:
			self.Dim = Dim

		self.dt = float(dt)
		self.N, self.dims = array(N[:Dim], int32), array(dims[:Dim], float64)
		self.h = self.dims / self.N
		self.rho, self.mu = float(rho), float(mu)
		self.DeltaType = DeltaType

		self.ConvectionMethod = ConvectionMethod
		
		ScalarDim, VectorDim, MatrixDim = self.FieldDims()
			
		self.u = zeros(VectorDim,float64)	 	# The fluid velocity vector field

		self.MakeWorkVars()

	def ConstitutiveName(self):
		return 'Stokes'
		
	def DimName(self):
		if self.Dim == 2:
			return '2D'
		else:
			return '3D'
    
	def FieldDims(self):
		N = self.N
		
		if self.Dim == 2:
			ScalarDim = (N[0],N[1])
			VectorDim = (N[0],N[1],2)
			MatrixDim = (N[0],N[1],2,2)
		else:
			ScalarDim = (N[0],N[1],N[2])
			VectorDim = (N[0],N[1],N[2],3)
			MatrixDim = (N[0],N[1],N[2],3,3)

		return ScalarDim, VectorDim, MatrixDim
		
	def MakeWorkVars(self):
		"""Create working variables used in the fluid solver."""
		
		ScalarDim, VectorDim, MatrixDim = self.FieldDims()
		
		self.Lookup = None # Used to store a lookup table for the Spread/Interpolated Stokes Greens function
		self.Symbols = Symbols.InitSymbols(self.N, self.h)

		self.du = zeros(MatrixDim,float64)		# The fluid total derivative (difference) matrix field
		self.p = zeros(ScalarDim,float64)		# The fluid pressure scalar field
		self.c = zeros(VectorDim,float64)	 	# The current timestep's explicit terms, a vector field
		self.f = zeros(VectorDim,float64)		# A vector field to potentially hold a force field

		self.u_Hat = zeros(VectorDim,complex64)	# The fluid velocity transform
		self.p_Hat = zeros(ScalarDim,complex64)	# The fluid pressure transform
		self.c_Hat = zeros(VectorDim,complex64)	# The transform of the explicit terms

		self.tempvec_complex = zeros(VectorDim,complex64)	# A vector field for storing temporary values
		self.tempscl_complex = zeros(ScalarDim,complex64)	# A scalar field for storing temporary values

		self.tempvec_real = zeros(VectorDim,float64)	# A vector field for storing temporary values
		self.tempscl_real = zeros(ScalarDim,float64)	# A scalar field for storing temporary values

		self.Output_u = self.u.copy()	# The output of the fluid solver, not associated with the actual evolution of the fluid
		
	def GetGrid(self):
		"""Returns a grid storing the indexical coordinates of each point in the Eulerian grid."""
		
		Dim, N = self.Dim, self.N
		
		if Dim == 2:
			return mgrid[0:N[0]:1, 0:N[1]:1]
		elif Dim == 3:
			return mgrid[0:N[0]:1, 0:N[1]:1, 0:N[2]:1]
		else:
			raise Exception("Domain must be either 2D or 3D")		
		
	def GetCoordinates(self):
		"""Return a grid storing the spatial coordinates of each point in the Eulerian grid."""
		
		Dim, dims, N = self.Dim, self.dims, self.N
		
		grid = self.GetGrid()
		coordinates = zeros(grid.shape, float64)
		
		for i in range(Dim):
			coordinates[i] = grid[i] * dims[i] / float(N[i])
			
		return coordinates
		
	def GetLookup(self):
		"""Returns the lookup table associated with the fluid."""
		
		if self.Lookup == None:
			self.Lookup = GetLookupTable(self)

		return self.Lookup
	
	def FluidSolve(self, f):
		"""Calculate one timeless fluid solve, with input force field f.
		The current fluid values are ignored, as is convection.
		The resulting fluid velocity is stored in Output_u."""
				
		self.c = self.FluidSolveExplicitTerms(f)

		FieldFFT(self.c, self.c_Hat)
	
		self.p_Hat = self.SolveFor_p_Hat()
				
		self.u_Hat = self.SolveFor_u_Hat()

		FieldIFFT(self.u_Hat, self.Output_u)
	
	def UpdateFluid(self, f, Update_du = True, CalcPressure = False, KeepUpdate = True):
		"""Advance the fluid one timestep via a semi-implicit scheme.
		
		f is a force field applied to the fluid at the current time and must be a vector field.
		if KeepUpdate is False then the updated fluid velocity is stored in self.Output_u, not self.u"""
	
		if Update_du:
			du = self.VectorGradient(self.u, self.du)
			
		self.c = self.ExplicitTerms(f)

		FieldFFT(self.c, self.c_Hat)
	
		self.p_Hat = self.SolveFor_p_Hat()
		
		self.u_Hat = self.SolveFor_u_Hat()

		if KeepUpdate:
			FieldIFFT(self.u_Hat, self.u)
			self.u = self.u.copy()
		else:
			FieldIFFT(self.u_Hat, self.Output_u)
			self.Output_u = self.Output_u.copy()
		
		if CalcPressure:
			self.p = fft.ifftn(self.p_Hat)

	def SolveFor_p_Hat(self):
		"""Solve for the pressure p in Fourier space via the formula
		
		p_Hat = D_Hat * c_Hat / (dt * Lambda_Hat)     <--- For Navier-Stokes
		p_Hat = D_Hat * c_Hat / (Lambda_Hat)          <--- For Stokes
		
		c_Hat must be the transform of explicit terms
		Lambda_Hat, D_Hat are the Fourier symbols of a discrete Laplacian and total difference operator respectively"""
		
		dt, Lambda_Hat, D_Hat, c_Hat, temp = self.dt, self.Symbols.WLambda_Hat, self.Symbols.D_Hat, self.c_Hat, self.tempscl_complex

		if self.rho == 0:
			p_Hat = FieldDot(D_Hat, c_Hat, temp, Dim = self.Dim) / Lambda_Hat			
		else:
			p_Hat = FieldDot(D_Hat, c_Hat, temp, Dim = self.Dim) / (dt * Lambda_Hat)

			
		if self.Dim == 2:
			p_Hat[0,0] = 0.
		else:
			p_Hat[0,0,0] = 0.
		
		return p_Hat

	def SolveFor_u_Hat(self):
		"""Solve for the velocity u in Fourier space via the formula
		
		u_Hat = (-dt * D_Hat * p_Hat + c_Hat) / (rho - dt * mu * Lambda_Hat) <--- For Navier-Stokes
		u_Hat = (-D_Hat * p_Hat + c_Hat) / (-mu * Lambda_Hat)                 <--- For Stokes
		
		rho, mu are fluid constants, the density and viscosity respectively
		c_Hat must be the transform of explicit terms
		p_Hat must be the transform of the associated pressure p
		Lambda_Hat, D_Hat are the Fourier symbols of a discrete Laplacian and total difference operator respectively"""

		dt, Lambda_Hat, D_Hat, p_Hat, c_Hat, temp, rho, mu = self.dt, self.Symbols.SLambda_Hat, self.Symbols.D_Hat, self.p_Hat, self.c_Hat, self.tempvec_complex, self.rho, self.mu

		if rho == 0:
			return FieldDot((-FieldDot(D_Hat, p_Hat, temp, Dim = self.Dim) + c_Hat), 1. / (-mu * Lambda_Hat), temp, Dim = self.Dim)
		else:	
			return FieldDot((-dt * FieldDot(D_Hat, p_Hat, temp, Dim = self.Dim) + c_Hat), 1. / (rho - dt * mu * Lambda_Hat), temp, Dim = self.Dim)

	def ExplicitTerms(self, f):
		"""Calculate explicit terms to be used in updating the fluid.
		
		Combines a force field f with convection rho * du * u and the current velocity u via the formula
		
		rho * u + dt * (f - rho * du * u)    <--- For Navier-Stokes
		f                                    <--- For Stokes
		
		du is a 3x3 array storing the central differences of u
		temp should be a vector field, it is used to hold temporary values"""

		dt, u, du, temp, rho = self.dt, self.u, self.du, self.tempvec_complex, self.rho
		
		if rho == 0:
			return f
		else:
			if self.ConvectionMethod == ConvectionMethods.ENO2:
				temp = self.Convect_Eno2(self.u)
			elif self.ConvectionMethod == ConvectionMethods.ENO3:
				temp = self.Convect_Eno3(self.u)
			elif self.ConvectionMethod == ConvectionMethods.Centered:
				FieldDot(du, u, temp, Dim = self.Dim)
			else:
				raise Exception("Unknown convection method specified.")
				
			temp *= -rho
			temp += f;
			temp *= dt
			val = rho * u
			val += temp
			
			return val

	def FluidSolveExplicitTerms(self, f):
		"""Calculates explicit terms to be used in a fluid solve.
		
		There is no convection and the current fluid velocity field is ignored.
		The resulting field is given by
		
		dt * f    <--- For Navier-Stokes
		f         <--- For Stoke's
		"""
		
		if self.rho == 0:
			return f
		else:
			dt = self.dt
		
			return dt * f

	def VectorGradient(self, f, df=None):
		"""Calculate the gradient of f and store it in df, if given."""
		
		if df == None:
			df = 0. * self.du
			
		Derivatives.TotalDifference(self.N, self.h, f, df, Derivatives.DifferenceTypes.CentralPeriodic)
		
		return df
	
	def DeformationTensor(self):
		"""Calculate the symmetric deformation tensor of the fluid."""
		
		self.du = self.VectorGradient(self.u, self.du)
		
		d = self.du * 0
		
		for i in range(self.Dim):
			for j in range(self.Dim):
				d[...,i,j] = .5 * (self.du[...,i,j] + self.du[...,j,i])
				
		return d

	def Convect_Eno2(self, phi):
		if self.Dim == 2:
			return self.Convect_Eno_2D(phi, Eno2)
		else:
			raise Exception("Not implemented for 3D")
		
	def Convect_Eno3(self, phi):
		if self.Dim == 2:
			return self.Convect_Eno_2D(phi, Eno3)
		else:
			raise Exception("Not implemented for 3D")
		
	def Convect_Eno_2D(self, phi, eno):
		f = 0 * phi
		u, N, h, Dim = self.u, self.N, self.h, self.Dim
		
		for k in range(phi.shape[Dim]):
			for i in range(N[0]):
				temp = f[i,:,1] * 0					
				eno(N[1], h[1], u[i,:,1], phi[i,:,k], temp)
				f[i,:,k] += temp
			for i in range(N[1]):
				temp = f[:,i,0] * 0				
				eno(N[0], h[0], u[:,i,0], phi[:,i,k], temp)
				f[:,i,k] += temp
				
		return f
			
	def Vorticity(self):
		"""Return the vorticity of the velocity field u."""
		
		du = self.VectorGradient(self.u, self.du)
		
		if self.Dim == 2:
			return du[...,1,0] - du[...,0,1]
		else:
			w = 0 * self.u
			w[...,0] = du[...,2,1] - du[...,1,2]
			w[...,1] = du[...,0,2] - du[...,2,0]
			w[...,2] = du[...,1,0] - du[...,0,1]
			
			return w
		
	def Vorticity2dSlice(self):
		"""Return a 2D slice of the vorticity field, dotted with the normal of the slice."""
		
		du = self.VectorGradient(self.u, self.du)
		
		if self.Dim == 2:
			return du[...,1,0] - du[...,0,1]
		else:
			return du[...,self.N[2]/2, 1,0] - du[...,self.N[2]/2, 0,1]
		
	def Params_ToString(self):
		"""Construct a string uniquely identifying the fluid parameters.
		Used to create unique file names."""
		
		string = 'N%g_dt%g' % (self.N[0], self.dt)
		
		if self.rho != 1 or self.mu != 1:
			string += '_rho%g_mu%g' % (self.rho, self.mu)
		
		return string
		
		#string = 'N_' + str(self.N[0]) + '_dt_' + str(self.dt).replace('.', 'P')

		#if self.rho != 1 or self.mu != 1:
			#string += '_rho_' + str(self.rho) + '_mu_' + str(self.mu)

		#return string

	def Extent(self):
		"""Return the extents of the fluid domain, as a 2 * Dim vector."""
		
		dims = self.dims
		
		if self.Dim == 2:
			return [0, dims[0], 0, dims[1]]
		else:
			return [0, dims[0], 0, dims[1], 0, dims[2]]
		
	def PlotExtent2D(self):
		"""Return the extents of a 2D slice of the fluid domain."""
		
		dims = self.dims
		return [0, dims[0], 0, dims[1]]
	
	def SetExtent(self):
		"""Sets the extents of the fluid domain."""
		
		dims = self.dims
		
		if self.Dim == 2:
			m.xlim(0, dims[0])
			m.ylim(0, dims[1])
		else:
			m.xlim(0, dims[0])
			m.ylim(0, dims[1])
		
		
		
	def SetPlotExtent(self):
		"""Sets the extent of the current plot to match the fluid domain."""
		
		extent = self.Extent()
		
		m.xlim(extent[0], extent[1])
		m.ylim(extent[2], extent[3])
		

	def PlotField(self, f):
		"""Plot a field f over the domain of the fluid."""
		
		self.Plot2D(field = f)
		
	def Plot2D(self, plot = PlotStyle.magnitude, Colorbar=True, field=None, additional=None, mod=None, colorbar_kwargs={}):
		"""Plot the current fluid.
		If 'field' is specified then 'field' is plotted, otherwise a component or function of the fluid is plotted,
		acoording to the plot parameter."""
		
		if field == None:
			if plot == PlotStyle.pressure:
				f = self.p.real
			elif plot == PlotStyle.vorticity:
				f = self.Vorticity()
			else:
				f = self.u
		else:
			f = field
		
		if self.Dim == 3:
			# The domain is 3D so take a cross section
			f = f[:,:,self.N[2]/2,:]
			
			if plot == PlotStyle.vorticity:
				f = f[...,2]
		
		if field != None:
			val = f
		elif plot == PlotStyle.magnitude:
			val = FieldDot(f, f, Dim = self.Dim)**.5
		elif plot == PlotStyle.xvel:
			val = f[...,0]
		elif plot == PlotStyle.yvel:
			val = f[...,1]
		elif plot == PlotStyle.zvel:
			val = f[...,2]
		elif plot == PlotStyle.pressure:
			val = f
		elif plot == PlotStyle.vorticity:
			val = f
		elif plot == PlotStyle.vel:
			self.PlotVectorField(f)
			return
		
		if mod != None: val = mod(val)
		m.imshow(val.T, interpolation='bilinear', extent=self.Extent(), origin='lower')
		#m.contour(val.T, interpolation='bilinear', extent=self.Extent(), origin='lower')

		if Colorbar: m.colorbar(**colorbar_kwargs)
		if additional != None: additional()
		
	def PlotVectorField(self, f):
		coords = self.GetCoordinates()
				
		N = self.N
		#PlotN = min(32, N[0])
		PlotN = min(24, N[0])
		
		coords = coords[:,::N[0]/PlotN,::N[1]/PlotN]
		f = f[::N[0]/PlotN,::N[1]/PlotN]
		
		#m.quiver(coords[0], coords[1], f[...,0], f[...,1], pivot='mid', color='r')
		m.quiver(coords[0], coords[1], f[...,0], f[...,1], pivot='mid')
		
	def __getstate__(self):
		d = dict(self.__dict__)
		del d['Lookup'], d['du'], d['f'], d['Symbols'], d['p'], d['c'], d['u_Hat'], d['p_Hat'], d['c_Hat'], d['tempvec_complex'], d['tempscl_complex'], d['Output_u'], d['tempscl_real'], d['tempvec_real']
		

		return d
	
	def __setstate__(self, dict):
		self.__dict__ = dict
	
		self.MakeWorkVars()
	