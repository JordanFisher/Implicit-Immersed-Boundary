from numpy import dot, array
from IBUtility import EnsureNumpyArray, cross
from Tethered import Tethered

class SquarePlate(Tethered):
	"""A 2D plate"""
	
	@staticmethod
	def N_to_dt(N, s):
		"""Return the largest stable timestep for an explicit simulation."""
		
		h = 1. / N
		return 10 * h / s**.5
		
	@staticmethod
	def N_to_Nb_x(N):
		"""Return the preferred value of Nb_x for a given value of N.
		
		Nb_x is the number of fiber points along one edge of the plate."""
		
		return int(N / 2 * 2.**.5)
	
	def __init__(self, Nb_x=0, s=0, p1=[.25,25,.5], p2=[.75,.25,.5], p3=[.75,.75,.5], center=None, size=None, h=None):
		"""Nb_x gives the number of fiber points along one side of the plate.
		
		s is an elasticity constant.
		p1, p2, p3 are sequential corners of the square.
		
		Alternatively, center specifies the center of the plate and size[0] and size[1] specify the widths in the x and y directions.
		
		h is a height function used to the perturb the height of the plate."""

		Tethered.__init__(self, (Nb_x + 1)**2, Dim = 3, s = s)

		if center != None:
			p1 = center + array([-size[0], -size[1], 0]) / 2
			p2 = center + array([size[0], -size[1], 0]) / 2
			p3 = center + array([size[0], size[1], 0]) / 2
		
		self.MakeSquare(Nb_x + 1, p1, p2, p3, h)

	def MakeSquare (self, N, p1, p2, p3, h):
		"""Make the fiber into a square sheet.
		
		p1, p2, p3 are sequential corners of the square.
		h is a height function used to the perturb the height of the plate."""
		
		if h == None: h = lambda x, y: 0
		
		p1, p2, p3 = EnsureNumpyArray(p1), EnsureNumpyArray(p2), EnsureNumpyArray(p3)

		self.Tether = self.X.copy()

		x, y = p2 - p1, p3 - p2
		z = cross(x, y)
		z /= dot(z,z)**.5

		hb = 1. / (N - 1)
		d1, d2 = x * hb, y * hb
		for j in range(N):
			for k in range(N):
				self.X[j + k * N,:] = p1 + j * d1 + k * d2 + z * h(j * hb, k * hb)
				self.Tether[j + k * N,:] = p1 + j * d1 + k * d2