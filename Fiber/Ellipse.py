from Fiber import Fiber

from numpy import cos, sin, pi
from Derivatives import Difference_SecondPeriodic1D
from IBUtility import EnsureNumpyArray

class Ellipse(Fiber):
	"""A 1D fiber in the shape of an ellipse.
	
	The force distribution is given by the second derivative times a stiffness constant."""
	
	def __init__(self, Nb, s, c, p1, p2, Dim = 2):
		"""Initialize an ellipsoidal fiber.
				
		Nb is the number of fiber points.
		s is a stiffness constant.
		
		c is the center.
		p1, p2 are the endpoints of the semimajor and minor axis."""
		
		Fiber.__init__(self, Nb, Dim)
		
		self.s = s
		
		self.MakeEllipse(c, p1, p2)

	def CalcForceDistribution(self):
		"""Calculate the force distribution along the fiber"""
		
		Difference_SecondPeriodic1D(self.Nb, self.hb, self.X, self.F)
		self.F *= self.s

		return self.F
		
	def CalcLinearForce(self):
		"""Calculate the linearized force."""
		
		return self.CalcForceDistribution()
		
	def MakeEllipse(self, c, p1, p2):
		"""Make the fiber into an ellipse
		c is the center
		p1, p2 are the endpoints of the semimajor and minor axis"""

		c, p1, p2 = EnsureNumpyArray(c), EnsureNumpyArray(p1), EnsureNumpyArray(p2)

		d1, d2 = p1 - c, p2 - c
		da = 2 * pi / self.Nb
		for i in range(self.Nb):
			a = i * da
			self.X[i] = cos(a) * d1 + sin(a) * d2 + c
