from Fluid import Fluid, ConvectionMethods
import Derivatives

from FieldUtility import FieldDot

import pylab as m
from numpy import array, zeros, float64, maximum, minimum
	
class ConstitutiveEqns:
	"""Different choices for the constitutive equation."""
	
	OB_TeranFauciShelley, OB_Chrispell, FENE_P = range(3)
	Name = { OB_TeranFauciShelley : 'OB_TFS', OB_Chrispell : 'OB_Chrispell', FENE_P : 'FENE_P' }
	
class UpdateMethods:
	ForwardEuler, RungeKutta2, RungeKutta3_TVD = range(3)
	
class SymmetricHelper2x2:
	_sym = array([[0,1], [1,2]]);

	_00 = _xx = 0
	_01 = _xy = 1
	_10 = _yx = 1
	_11 = _yy = 2
	
	UniqueEntries = 3

	@staticmethod
	def identity():
		"""Return the identity matrix."""
		
		return array([ 1.,  0.,  1.])

	@classmethod
	def FieldDot(cls, S, x, row, y):
		"""At every point of the domain calculate the dot product of a row of S with x, storing the result in y.
		S is a symmetric matrix field stored in an array of shape NxNxNx6"""
		
		_sym = cls._sym

		y *= 0

		for i in range(2):
			y[...] += S[...,_sym[row,i]] * x[...,i]

		return y

	@classmethod
	def Trace(cls, S):
		"""Take a tensor field S and return the trace, a scalar field."""
		
		sym = cls
		
		return S[...,sym._00] + S[...,sym._11]
	
	ComponentNames = ['xx', 'xy', 'yy']
		
class SymmetricHelper3x3:
	_sym = array([[0,1,2], [1,3,4], [2,4,5]]);

	_00 = _xx = 0
	_01 = _xy = 1
	_02 = _xz = 2
	_10 = _yx = 1
	_11 = _yy = 3
	_12 = _yz = 4
	_20 = _zx = 2
	_21 = _zy = 4
	_22 = _zz = 5

	UniqueEntries = 6
	
	@staticmethod
	def identity():
		"""Return the identity matrix."""
		
		return array([ 1.,  0.,  0.,  1.,  0.,  1.])
	
	@classmethod
	def FieldDot(cls, S, x, row, y):
		"""At every point of the domain calculate the dot product of a row of S with x, storing the result in y.
		S is a symmetric matrix field stored in an array of shape NxNxNx6"""
		
		_sym = cls._sym

		y *= 0

		for i in range(3):
			y[...] += S[...,_sym[row,i]] * x[...,i]

		return y

	@classmethod
	def Trace(cls, S):
		"""Take a tensor field S and return the trace, a scalar field."""
		
		sym = cls
		
		return S[...,sym._00] + S[...,sym._11] + S[...,sym._22]

	ComponentNames = ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']
	
class ViscoFluid(Fluid):
	"""Holds all the properties and data associated to a boxed fluid governed by the Oldroy-B model"""

	def __init__ (self, dt, N, dims, W = 1., B = .5, L = 100., rho = 1., mu = 1., DeltaType = 0, Dim = None,
	              ConstitutiveEqn = ConstitutiveEqns.OB_TeranFauciShelley,
	              ConvectionMethod = ConvectionMethods.Upwind):
		"""dims should be a 2-vector or 3-vector holding the dimensions of the fluid box.
		N should likewise be an integer 2-vector or 3-vector which describes how fine the grid will be in each dimension.
		Dim is the number of spatial dimensions (2 or 3) of the fluid domain."""

		Fluid.__init__(self, dt, N, dims, rho, mu, DeltaType = DeltaType, Dim = Dim)
		
		self.ConstitutiveEqn = ConstitutiveEqn
		self.ConvectionMethod = ConvectionMethod
		
		Dim = self.Dim
		
		# Get the relevant symmetric matrix helper class
		if Dim == 2:
			self.MySymmetricHelper = SymmetricHelper2x2
		elif Dim == 3:
			self.MySymmetricHelper = SymmetricHelper3x3
		else:
			raise Exception("Dim, the dimension of the fluid domain, must be either 2 or 3.")
		
		# The various array dimensions for scalar/vector/matrix fields, depending on the spatial dimension of the fluid domain.
		if Dim == 2:
			ScalarDim = (N[0],N[1])
			VectorDim = (N[0],N[1],2)
			SymMatrixDim = (N[0],N[1],3)
			dSDim = (N[0],N[1],3,2)
		else:
			ScalarDim = (N[0],N[1],N[2])
			VectorDim = (N[0],N[1],N[2],3)
			SymMatrixDim = (N[0],N[1],N[2],6)
			dSDim = (N[0],N[1],N[2],6,3)
		
		self.W = float(W)  # The Weissenberg number
		self.B = float(B)  # beta
		self.L = float(L)  # Maximum extension length for FENE
		self.S = zeros(SymMatrixDim,float64) # The polymeric stress (or conformation tensor), with 6 components at each material point (a 3x3 symmetric matix)
		self.CurrentStress = None # The actual polymeric stress, calculated from the conformation tensor
		
		# Initialize the stress (or conformation tensor)
		if self.ConstitutiveEqn == ConstitutiveEqns.OB_Chrispell or self.ConstitutiveEqn == ConstitutiveEqns.OB_TeranFauciShelley:
			self.S[...] =  self.MySymmetricHelper.identity() # Initial value is the 2x2 (3x3) identity matrix
		else:
			self.S[...] =  self.MySymmetricHelper.identity() # Initial value is the 2x2 (3x3) identity matrix
			#self.S *= 0 # Initial value is the zero matrix
		
		self.dS = zeros(dSDim,float64) # The total derivative of S

		self.g = zeros(SymMatrixDim,float64) # Used to update S
		
		self.ArtificialDiffusion = self.h[0]**2

	def TensorGradient(self, f, df=None):
		"""Calculate the gradient of a symmetric tensor field f and store it in df, if given."""
		
		if df == None:
			df = 0. * self.dS
			
		Derivatives.TotalDifference(self.N, self.h, f, df, Derivatives.DifferenceTypes.CentralPeriodic)
		
		return df
	
	def LeftRightTensorGradient(self, f):
		"""Calculate the left and right difference gradient of a symmetric tensor field f."""
		
		df_Left, df_Right = 0. * self.dS, 0. * self.dS
			
		Left = Derivatives.TotalDifference(self.N, self.h, f, df_Left, Derivatives.DifferenceTypes.LeftPeriodic)
		Right = Derivatives.TotalDifference(self.N, self.h, f, df_Right, Derivatives.DifferenceTypes.RightPeriodic)		
		
		return df_Left, df_Right

	def StressTrace(self):
		"""Get the trace of S."""
		
		sym = self.MySymmetricHelper
		
		trace = sym.Trace(self.S)

		return trace
	
	def PeterlinFunc(self, S):
		"""Calculate the Peterlin function, used in the FENE-P model.
		l is the maximum length of the FENE dumbbell,
		S is the conformation tensor field.
		
		A scalar field is returned, where each point in the field is given by the formula
		l / (l - trace(S(x)))"""
		
		sym = self.MySymmetricHelper

		L = self.L
		
		trace = sym.Trace(S)
		
		return L / (L - trace)
	
	def Calc_g2D(self, S, dS=None, ConvectionMethod=None):
		"""g = -(u * \nabla)S +
			   \nabla u * S + S * \nabla u^T +
			   1/W * (I - S)"""
			   
		if dS == None:
			dS = self.TensorGradient(S, dS)

		if ConvectionMethod == None:
			ConvectionMethod = self.ConvectionMethod
		
		g = self.g
		u = self.u
		du = self.du
		W = self.W
		B = self.B
		Dim = self.Dim

		sym = self.MySymmetricHelper

		sDot = sym.FieldDot

		# Calculate the convection
		if ConvectionMethod == ConvectionMethods.ENO2:
			convection = self.Convect_Eno2(S)
			g[...] = -convection
		elif ConvectionMethod == ConvectionMethods.ENO3:
			convection = self.Convect_Eno3(S)
			g[...] = -convection			
		elif ConvectionMethod == ConvectionMethods.Upwind:
			dS_Left, dS_Right = self.LeftRightTensorGradient(S)

			uMax = maximum(u, 0)
			uMin = minimum(u, 0)
			
			g[...,sym._00] = -dS_Left[...,sym._00,0] * uMax[...,0] - dS_Left[...,sym._00,1] * uMax[...,1]
			g[...,sym._01] = -dS_Left[...,sym._01,0] * uMax[...,0] - dS_Left[...,sym._01,1] * uMax[...,1]
			g[...,sym._11] = -dS_Left[...,sym._11,0] * uMax[...,0] - dS_Left[...,sym._11,1] * uMax[...,1]

			g[...,sym._00] += -dS_Right[...,sym._00,0] * uMin[...,0] - dS_Right[...,sym._00,1] * uMin[...,1]
			g[...,sym._01] += -dS_Right[...,sym._01,0] * uMin[...,0] - dS_Right[...,sym._01,1] * uMin[...,1]
			g[...,sym._11] += -dS_Right[...,sym._11,0] * uMin[...,0] - dS_Right[...,sym._11,1] * uMin[...,1]
		elif ConvectionMethod == ConvectionMethods.Centered:
			g[...,sym._00] = -dS[...,sym._00,0] * u[...,0] - dS[...,sym._00,1] * u[...,1]
			g[...,sym._01] = -dS[...,sym._01,0] * u[...,0] - dS[...,sym._01,1] * u[...,1]
			g[...,sym._11] = -dS[...,sym._11,0] * u[...,0] - dS[...,sym._11,1] * u[...,1]
		else:
			raise Exception("Unknown convection method specified.")
		
		# Additional component from upper convected derivative (frame stretching)
		g[...,sym._00] += 2 * (du[...,0,0] * S[...,sym._00] + du[...,0,1] * S[...,sym._01])
		g[...,sym._01] +=      du[...,0,1] * S[...,sym._11] + du[...,1,0] * S[...,sym._00]
		g[...,sym._11] += 2 * (du[...,1,0] * S[...,sym._01] + du[...,1,1] * S[...,sym._11])
	
		#Teran-Fauci-Shelley's equation
		if self.ConstitutiveEqn == ConstitutiveEqns.OB_TeranFauciShelley:
			g[...,sym._00] += (1. / W) * (1 - S[...,sym._00])
			g[...,sym._01] += (1. / W) * (0 - S[...,sym._01])
			g[...,sym._11] += (1. / W) * (1 - S[...,sym._11])
		# Chrispell's equation
		elif self.ConstitutiveEqn == ConstitutiveEqns.OB_Chrispell:
			d = self.DeformationTensor()		
	
			g[...,sym._00] += (1. / W) * (2 * B * d[...,0,0] - S[...,sym._00])
			g[...,sym._01] += (1. / W) * (2 * B * d[...,0,1] - S[...,sym._01])
			g[...,sym._11] += (1. / W) * (2 * B * d[...,1,1] - S[...,sym._11])
		# FENE-P
		elif self.ConstitutiveEqn == ConstitutiveEqns.FENE_P:
			f = self.PeterlinFunc(S)
			
			g[...,sym._00] += (1. / W) * (1 - f * S[...,sym._00])
			#z=2 * (du[...,0,0] * S[...,sym._00] + du[...,0,1] * S[...,sym._01])
			#w=(1. / W) * (1 - f * S[...,sym._00])
			#n=self.N[1]/2
			#m=S[:,n,0].argmax()
			#z=z[m,n]
			#w=w[m,n]
			#s=S[m,n,0]
			#print s, z, w
			g[...,sym._01] += (1. / W) * (0 - f * S[...,sym._01])
			g[...,sym._11] += (1. / W) * (1 - f * S[...,sym._11])
		else:
			raise Exception("Unknown Oldroy-B model.")
		
		# Add artificial dissipation
		g += self.ArtificialDiffusion * Derivatives.VectorLaplacian(self.N, self.h, S)

		return g		

	def Calc_g3D(self, S, dS=None, ConvectionMethod=None):
		"""g = -(u * \nabla)S +
			   \nabla u * S + S * \nabla u^T +
			   1/W * (I - S)"""
			   
		if dS == None:
			dS = self.TensorGradient(S, dS)

		if ConvectionMethod == None:
			ConvectionMethod = self.ConvectionMethod
		
		g = self.g
		u = self.u
		du = self.du
		W = self.W
		B = self.B
		Dim = self.Dim

		sym = self.MySymmetricHelper

		sDot = sym.FieldDot

		# Calculate the convection
		if ConvectionMethod == ConvectionMethods.Upwind:
			dS_Left, dS_Right = self.LeftRightTensorGradient(S)

			uMax = maximum(u, 0)
			uMin = minimum(u, 0)

			for i in range(sym.UniqueEntries):
				g[...,i] =  -dS_Left[...,i,0] * uMax[...,0] - dS_Left[...,i,1] * uMax[...,1] - dS_Left[...,i,2] * uMax[...,2]
				g[...,i] += -dS_Right[...,i,0] * uMin[...,0] - dS_Right[...,i,1] * uMin[...,1] - dS_Right[...,i,2] * uMin[...,2]
		elif ConvectionMethod == ConvectionMethods.Centered:
			for i in range(sym.UniqueEntries):
				g[...,i] =  -dS[...,i,0] * u[...,0] - dS[...,i,1] * u[...,1] - dS[...,i,2] * u[...,2]
		else:
			raise Exception("Unknown convection method specified.")
	
		# Additional component from upper convected derivative (frame stretching)
		g[...,sym._00] += 2 * sDot(S[...], du[...,0,:], 0, self.tempscl_real)
		g[...,sym._11] += 2 * sDot(S[...], du[...,1,:], 1, self.tempscl_real)
		g[...,sym._22] += 2 * sDot(S[...], du[...,2,:], 2, self.tempscl_real)

		g[...,sym._01] += sDot(S[...], du[...,0,:], 1, self.tempscl_real) + sDot(S[...], du[...,1,:], 0, self.tempvec_real[...,0])
		g[...,sym._02] += sDot(S[...], du[...,0,:], 2, self.tempscl_real) + sDot(S[...], du[...,2,:], 0, self.tempvec_real[...,0])
		g[...,sym._12] += sDot(S[...], du[...,1,:], 2, self.tempscl_real) + sDot(S[...], du[...,2,:], 1, self.tempvec_real[...,0])
	
		#Teran-Fauci-Shelley's equation
		if self.ConstitutiveEqn == ConstitutiveEqns.OB_TeranFauciShelley:
			I = sym.identity()
			for i in range(sym.UniqueEntries):
				g[...,i] += (1. / W) * (I[i] - S[...,i])
		# FENE-P
		elif self.ConstitutiveEqn == ConstitutiveEqns.FENE_P:
			f = self.PeterlinFunc(S)
			
			I = sym.identity()
			for i in range(sym.UniqueEntries):
				g[...,i] += (1. / W) * (I[i] - f * S[...,i])
		else:
			raise Exception("Unknown Consitutive equation.")


		# Add artificial dissipation
		g += self.ArtificialDiffusion * Derivatives.VectorLaplacian(self.N, self.h, S)

		return g
		
	
	def ImplicitFENE_PStep(self):
		"""Experimental FENE-P implicit Euler step, make sure to modify Calc_g."""
		
		sym = self.MySymmetricHelper
		W, L = self.W, self.L
		
		def E(V):
			A = V + dt * self.Calc_g(V)
			trA = sym.Trace(A)

			# Solve the non-linear system (which amounts to a quadratic)
			a = 1.
			b = -(trA + (1. + dt / W) * L)
			c = trA * L
			discriminant = b**2 - 4 * a * c
			trS = .5 * (-b - discriminant**.5)
			
			check = trA - (dt * L / W) * trS / (L - trS) - trS
			print "check", abs(check).max()
			print "desired trace", trS.max(), trS.min()
			oldA = A * 1.
			
			# S = A - factor * S
			factor = (dt / W * L) / (L - trS)
			print factor.max(), factor.min()
			A[...,sym._00] = A[...,sym._00] / (1 + factor)
			A[...,sym._01] = A[...,sym._01]
			A[...,sym._11] = A[...,sym._11] / (1 + factor)
			
			newtr = sym.Trace(A)
			print "new trace", newtr.max(), newtr.min()
			check = oldA[...,sym._00] - (dt * L / W) * A[...,sym._00] / (L - newtr) - A[...,sym._00]
			print "final check", abs(check).max()
			check = oldA[...,sym._11] - (dt * L / W) * A[...,sym._11] / (L - newtr) - A[...,sym._11]
			print "final check", abs(check).max()

			
			return A
		
		return E
		
	def UpdateS(self, Steps = 1, Method = UpdateMethods.RungeKutta2):
		"""Update the stress S.
		
		The timestep taken is the fluid timstep divided by Steps. This can increase stability."""
		
		dt = self.dt / Steps

		# Define the forward Euler function
		# Various Runge-Kutta methods can be written in terms of this function
		if self.Dim == 2:
			def E(V):
				return V + dt * self.Calc_g2D(V)
		else:
			def E(V):
				return V + dt * self.Calc_g3D(V)
			
		# Experimental FENE-P implicit Euler step
		#E = self.ImplicitFENE_PStep()
		#Method = UpdateMethods.ForwardEuler
			
		for i in range(Steps):
			# Forward Euler
			if Method == UpdateMethods.ForwardEuler:
				self.S = E(self.S)
			
			# Runge-Kutta 2
			elif Method == UpdateMethods.RungeKutta2:
				S1 = E(self.S)		   
				S2 = E(S1)
				
				self.S = (self.S + S2) / 2.

			# Runge-Kutta 3, TVD
			elif Method == UpdateMethods.RungeKutta3_TVD:
				self.S = 1./3 * S + 2./3 * E(3./4 * S + 1./4 * E(E(self.S)))
	
		
	def GetTau(self):
		"Gets the stress tensor tau."""

		sym = self.MySymmetricHelper

		# Calculate the stress tensor
		if self.ConstitutiveEqn == ConstitutiveEqns.OB_TeranFauciShelley:
			# \beta S
			tau = self.S
		elif self.ConstitutiveEqn == ConstitutiveEqns.OB_Chrispell:
			# S
			tau = self.S
		elif self.ConstitutiveEqn == ConstitutiveEqns.FENE_P:
			# (1/W) * (fS - I)
			f = self.PeterlinFunc(self.S)
			
			tau = 0 * self.S
			tau[...,sym._00] = (1. / self.W) * (f * self.S[...,sym._00] - 1.)
			tau[...,sym._01] = (1. / self.W) * (f * self.S[...,sym._01] - 0.)
			tau[...,sym._11] = (1. / self.W) * (f * self.S[...,sym._11] - 1.)
			
			if self.Dim == 3:
				tau[...,sym._02] = (1. / self.W) * (f * self.S[...,sym._02] - 0.)
				tau[...,sym._12] = (1. / self.W) * (f * self.S[...,sym._12] - 0.)
				tau[...,sym._22] = (1. / self.W) * (f * self.S[...,sym._22] - 1.)
		else:
			raise Exception("Unknown constitutive equation.")		

		return tau
		
	def PolymericStress(self):
		"""Calculate the polymer's contribution to the stress."""
		
		sym = self.MySymmetricHelper
		
		# Get tau, the stress tensor
		tau = self.GetTau()

		# Calculate the gradient of tau
		dtau = self.TensorGradient(tau)		
		
		# Calculate and return the divergence of tau
		if self.Dim == 2:
			self.tempvec_real[...,0] = dtau[...,sym._00,0] + dtau[...,sym._01,1]
			self.tempvec_real[...,1] = dtau[...,sym._10,0] + dtau[...,sym._11,1]
		else:
			self.tempvec_real[...,0] = dtau[...,sym._00,0] + dtau[...,sym._01,1] + dtau[...,sym._02,2]
			self.tempvec_real[...,1] = dtau[...,sym._10,0] + dtau[...,sym._11,1] + dtau[...,sym._12,2]
			self.tempvec_real[...,2] = dtau[...,sym._20,0] + dtau[...,sym._21,1] + dtau[...,sym._22,2]		

		# Multiplicative factor
		if self.ConstitutiveEqn == ConstitutiveEqns.OB_TeranFauciShelley:
			# \beta
			self.tempvec_real *= self.B
		elif self.ConstitutiveEqn == ConstitutiveEqns.OB_Chrispell:
			# 1
			self.tempvec_real *= 1
		elif self.ConstitutiveEqn == ConstitutiveEqns.FENE_P:
			# \beta
			self.tempvec_real *= self.B
		else:
			raise Exception("Unknown constitutive equation.")		

		return self.tempvec_real
				
	def UpdateFluid(self, f, Update_du = True, CalcPressure = False, UpdateStress = True, KeepUpdate = True):
		"""Advance the fluid one timestep via a semi-implicit scheme.
		
		f is a force field applied to the fluid at the current time and must be a vector field."""

		# Calculate polymeric stress contribution
		self.CurrentStress = self.PolymericStress()
		f += self.CurrentStress
		
		# Update the fluid
		Fluid.UpdateFluid(self, f, Update_du, CalcPressure, KeepUpdate)

		# Update the fluid stress
		if UpdateStress and KeepUpdate:
			self.UpdateS()


	def PlotComponent(self, val, name, subplot, Colorbar, f, mod):
		m.imshow(mod(val).T, interpolation='bilinear', extent=self.PlotExtent2D(), origin='lower')

		self.SetPlotExtent()
		
		# Suppress labels on the axis tips
		for tick in subplot.xaxis.get_ticklabels():
			tick.set_visible(False)
		for tick in subplot.yaxis.get_ticklabels():
			tick.set_visible(False)
		
		# Add a subtitle
		m.title(name)
		
		if Colorbar: m.colorbar(format='%0.5f')
		if f != None: f()

	def StressPlot(self, f=None, mod=None, Colorbar = True):
		"""Plot the components of the stress tensor.
		
		If the domain is 3D a 2D cross-section is plotted (z = N/2) and the xx/xy/yy/yz components are plotted.
		Otherwise the xx/xy/yx/yy components are plotted.
		
		f is an optional function that performs some additional action after plotting each subplot.
		mod is an optional function that takes a field as an argument and does some processing before plotting."""
		
		sym = self.MySymmetricHelper
		
		if mod == None:
			mod = lambda u: u

		tau = self.GetTau()
			
		if self.Dim == 2:
			Components = [sym._xx, sym._xy, sym._yx, sym._yy]
		else:
			tau = tau[:,:,self.N[2]/2]
			Components = [sym._xx, sym._xy, sym._yx, sym._yy]
	
		subplot = m.subplot(221)
		self.PlotComponent(tau[:,:,Components[0]], sym.ComponentNames[Components[0]] + ' component', subplot, Colorbar, f, mod)
		
		subplot = m.subplot(222)		
		vorticity = self.Vorticity2dSlice()
		self.PlotComponent(vorticity, 'Vorticity', subplot, Colorbar, f, mod)
		#self.PlotComponent(tau[:,:,Components[1]], sym.ComponentNames[Components[1]] + ' component', subplot, Colorbar, f, mod)
		
		subplot = m.subplot(223)
		self.PlotComponent(tau[:,:,Components[2]], sym.ComponentNames[Components[2]] + ' component', subplot, Colorbar, f, mod)
		
		subplot = m.subplot(224)
		self.PlotComponent(tau[:,:,Components[3]], sym.ComponentNames[Components[3]] + ' component', subplot, Colorbar, f, mod)
		
	def StressDebug(self, f=None, mod=None, Colorbar = True):
		sym = self.MySymmetricHelper
		
		if mod == None:
			mod = lambda u: u
		
		if self.Dim == 2:
			S = self.S
			Components = [sym._xx, sym._xy, sym._yx, sym._yy]
		else:
			S = self.S[:,:,self.N[2]/2]
			Components = [sym._xx, sym._xy, sym._yy, sym._yz]		
	
		subplot = m.subplot(221)
		self.PlotComponent(S[:,:,Components[0]], sym.ComponentNames[Components[0]] + ' component', subplot, Colorbar, f, mod)
		
		subplot = m.subplot(222)
		self.PlotComponent(self.Vorticity(), 'Vorticity', subplot, Colorbar, f, mod)
		
		subplot = m.subplot(223)
		self.PlotComponent(self.u[:,:,0], 'u', subplot, Colorbar, f, mod)
		
		subplot = m.subplot(224)
		self.PlotComponent(self.du[:,:,0,0], 'du/dx', subplot, Colorbar, f, mod)
		
	def SaveState(self, Saver):
		Fluid.SaveState(self, Saver)
		
		Saver.Save(self.S, 'S_' + str(Saver.CurStep))

	def __getstate__(self):
		d = Fluid.__getstate__(self)
		
		del d['g'], d['dS'], d['CurrentStress']
		
		return d