import pylab as m
from numpy import zeros, float64

from Delta import FiberToGrid, GridToFiber

import IB_c

class Fiber:
	"""An immersed fiber/boundary/structure.
	
	This is a generic class and does not define any specific geometry."""

	@classmethod
	def FromArray(cls, X):
		"""Create a fiber from an array of points X."""
		
		return cls(Nb = X.shape[0], Dim = X.shape[1])
		
	def __init__(self, Nb = 0, Dim = 3, UseNativePython = False):
		"""Set aside the needed memory to handle the fiber.
		
		Nb is the number of fiber points.
		Dim is the dimension of the underlying fluid domain.
		When UseNativePython is set to True, spreading and interpolation is perfomed in native python code."""

		self.Dim = Dim
		self.UseNativePython = UseNativePython
		
		self.X = zeros((Nb,Dim),float64)
		self.F = zeros((Nb,Dim),float64)
		self.U = zeros((Nb,Dim),float64)

		self.Nb = Nb
		self.hb = 1. / self.Nb

	"""Construct a matrix representation of the fluid interactions between all fiber points."""
	def ConstructFluidMatrix(self, fluid, M=None):
		if M == None:
			# Allocate memory to store the matrix
			Dim = fluid.Dim
			M = zeros((Dim * self.Nb, Dim * self.Nb), float64)
		else:
			M *= 0
			
		# Construct the matrix
		IB_c.ComputeMatrix(self.X, fluid.N, self.Nb, fluid.h, self.hb, fluid.GetLookup(), M, 100000)
	
		return M
		
	def Size(self):
		"""Return the size of the array holding X: Nb * Dim."""
		
		return self.Nb * self.Dim
		
	def SetTime(self, t):
		"""Update any time dependent properties of the fiber."""
		
		pass
		
	def SetLagged(self, X=None):
		"""Set the configuration used for spreading and interpolation."""
		
		if X == None: X = self.X
		
		self.LaggedX = X.copy()
		
	def CalcForceDistribution(self):
		"""Calculate the force distribution along the fiber."""

		raise Exception("Force distribution function not defined!")
		
	def CalcLinearForce(self):
		"""Calculate the linearized force."""
		
		raise Exception("Linear force function not defined!")
		
	def ToGrid(self, fluid, F=None, f=None, X=None):
		"""Spread a distribution defined on the fiber to the underlying grid.
		
		F is the defined distribution. The default value is the fiber's force distribution.
		f is the target Eulerian grid. The default target is fluid.f, a temporary vector field used to store fluid forces.
		X is the fiber configuration used. The default is the current configuration, self.X"""
		
		if F == None: F = self.F
		if f == None: f = fluid.f
		if X == None: X = self.X
		
		FiberToGrid (fluid.N, fluid.h, self.Nb, self.hb, X, F, f, fluid.DeltaType, self.UseNativePython)
			
	def FromGrid(self, fluid, u=None, U=None, X=None):
		"""Interpolate a vector field over the fluid domain to the fiber configuration.
		
		u is the vector field. The default value is the fluid's velocity field.
		U is the target Lagrangian vector field. The default target is self.U, a temporary vector field used to store velocity.
		X is the fiber configuration used. The default is the current configuration, self.X"""
		
		if u == None: u = fluid.u
		if U == None: U = self.U
		if X == None: X = self.X

		GridToFiber (fluid.N, fluid.h, self.Nb, self.hb, X, u, U, fluid.DeltaType, self.UseNativePython)
			
	def CheckForExplosion(self, Limit = 1):
		if abs(self.X).max() > Limit:
			print "Maximum fiber position norm is", abs(self.X).max()
			str = "Simulation has exploded: the immersed boundary has moved past a set limit."
				
			raise Exception(str)

	def Plot2D(self, **Scatter_args):
		"""Plot the current configuration."""
		
		X = self.X
		
		m.scatter(X[:,0], X[:,1], **Scatter_args)
	
		#m.xlim(0, 1)
		#m.ylim(0, 1)
		
	def PlotVel(self, **Scatter_args):
		"""Plot the current configuration."""
		
		X, U = self.X, self.U
		
		m.quiver(X[:,0], X[:,1], U[:,0], U[:,1], **Scatter_args)
	
		m.xlim(0, 1)
		m.ylim(0, 1)
		
	
	def SaveState(self, Saver):
		"""Save the current configuration of the fiber."""
		
		Saver.Save(self.X, 'X_' + str(Saver.CurStep))
