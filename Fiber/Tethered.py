from Fiber import Fiber

class Tethered(Fiber):
	"""A tethered immersed fiber/boundary/structure.
	
	This is a generic class and does not define any specific geometry."""
	
	def __init__(self, Nb = 0, Dim = 3, s = 1):
		"""Set aside the needed memory to handle the fiber."""
		
		Fiber.__init__(self, Nb, Dim)
		
		self.Tether = self.X.copy()
		self.s = s
		
	def CorrectStokes(self, t, fluid):
		"""Modify the fluid velocity field so that the fiber configuration yields a zero sum force."""

		self.SetTime(t + fluid.dt)
		
		self.F = (self.Tether - self.X)
		f = self.F.sum(axis=0) / self.Nb

		self.X += f

		fluid.u += f / fluid.dt
		
		
	def CalcForceDistribution(self):
		"""Calculate the force distribution along the fiber."""
		
		self.F = self.s * (self.Tether - self.X)
		
		return self.F
		
	def CalcLinearForce(self):
		"""Calculate the linearized force."""
		
		self.F = -self.s * self.X
		
		return self.F