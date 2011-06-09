from numpy import dot, array
from IBUtility import EnsureNumpyArray, cross
from FiberUtility import MakeTriangle, TrianglesToVEF
from Tethered import Tethered

class Sphere(Tethered):
	"""A hollow sphere, formed via tether points."""
	
	@staticmethod
	def N_to_dt(N, s):
		"""Return the largest stable timestep for an explicit simulation."""
		
		h = 1. / N
		return .05 * h / s**.5
		
	@staticmethod
	def N_to_Nb_x(N):
		"""Return the preferred value of Nb_x for a given value of N.
		
		Nb_x is the number of fiber points along the edge of one triangle of the generating icosohedron."""
		
		return int(N / 2 * 2.**.5)
	
	@staticmethod
	def MakeSphere(N, c, r):
		"""Create a discretized sphere centered at c with radius r.
		
		N determines how many points to use in the discretization."""
		
		c = EnsureNumpyArray(c)

		Triangles, Points = [], []

		# Make an icosohedron, where each face is further triangulated.
		MakeTriangle(N, c + [0,-r,0], c + [r,0,0], c + [0,0,r], Triangles = Triangles)
		MakeTriangle(N, c + [0,-r,0], c + [-r,0,0], c + [0,0,r], Triangles = Triangles)
		MakeTriangle(N, c + [0,r,0], c + [r,0,0], c + [0,0,r], Triangles = Triangles)
		MakeTriangle(N, c + [0,r,0], c + [-r,0,0], c + [0,0,r], Triangles = Triangles)
		MakeTriangle(N, c + [0,-r,0], c + [r,0,0], c + [0,0,-r], Triangles = Triangles)
		MakeTriangle(N, c + [0,-r,0], c + [-r,0,0], c + [0,0,-r], Triangles = Triangles)
		MakeTriangle(N, c + [0,r,0], c + [r,0,0], c + [0,0,-r], Triangles = Triangles)
		MakeTriangle(N, c + [0,r,0], c + [-r,0,0], c + [0,0,-r], Triangles = Triangles)

		# Convert the list of triangles into Verticies, Edges, and Triangles
		Points, Edges, Triangles = TrianglesToVEF(Triangles)

		Points = array(Points)

		# Project every point onto the sphere the icosohedron is inscribed in.
		# First shift the points so the center is the origin.
		dif = Points - c
		difsq = dif**2
		l = (difsq[:,0] + difsq[:,1] + difsq[:,2])**.5
		dif[:,0] *= r / l
		dif[:,1] *= r / l
		dif[:,2] *= r / l

		# Shift the points back to the original center.
		Points = dif + c

		return Points, array(Triangles), array(Edges)	
	
	def __init__(self, Nb_x=0, s=0, center=[.5,.5,.5], r=.2):
		"""Nb_x is the number of fiber points along the edge of one triangle of the generating icosohedron.
		
		s is an elasticity constant."""

		center = EnsureNumpyArray(center)
		
		# Make the sphere
		Points, Triangles, Edges = self.MakeSphere(Nb_x, center, r)
		
		# Make the fiber
		Tethered.__init__(self, Nb = Points.shape[0], Dim = 3, s = s)
		
		self.X = Points
		self.Tether = self.X.copy()
		
		self.s = s
		
		# Save the corresponding unit sphere
		self.UnitSphere = (self.X - center) / r
