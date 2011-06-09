from numpy import zeros, array, float64
from IBUtility import EnsureNumpyArray
import Fiber, Tethered

def TriangulateTriangle(Triangle):
	"""Take a single triangle (a 3-tuple of 3-vectors) and return an array of triangles which subdivide it."""
	
	Triangles = []
	a, b, c = Triangle
	ab, bc, ac = (a+b)/2, (b+c)/2, (a+c)/2
	return [(a, ab, ac), (b, ab, bc), (c, bc, ac), (ab, bc, ac)]

def Triangulate(Triangles):
	"""Take a list of triangles and return a new list where each triangle has been replaced by a subdivided triangle."""
	
	return [NewTriangle for Triangle in Triangles for NewTriangle in TriangulateTriangle(Triangle)]

def TrianglesToVEF(Triangles):
	"""Take a list of triangles and return a lists of Vertices, Edges, and Faces."""
	
	# Create a list of all the points in the list of triangles (including redundancies)
	Points = []
	for Triangle in Triangles:
		Points.append((Triangle[0][0], Triangle[0][1], Triangle[0][2]))
		Points.append((Triangle[1][0], Triangle[1][1], Triangle[1][2]))
		Points.append((Triangle[2][0], Triangle[2][1], Triangle[2][2]))		

	# Create a dictionary that takes a point and spits out a unique integer from 0 to the number of unique points minus 1
	PointDic = {}
	map(lambda x: PointDic.__setitem__(x, 0), Points)
	Points = PointDic.keys()
	i = 0
	for key in Points:
		PointDic[key] = i
		i += 1

	# Construct a new list of triangles where each vertex of a triangle is specified by its unique integer value in the point dictionary
	NewTriangles = [[PointDic[(Triangle[i][0], Triangle[i][1], Triangle[i][2])] for i in range(3)] for Triangle in Triangles]

	# Construct a list of edges where each vertex is specified by its unique integer value (this list of edges includes redundancies)
	Edges = []
	for Triangle in NewTriangles:
		for i in range(-1,2):
			if Triangle[i] < Triangle[i+1]:
				Edges.append( (Triangle[i], Triangle[i+1]) )
			else:
				Edges.append( (Triangle[i+1], Triangle[i]) )
				
	# Get rid of redundancies in the Edge list
	Edges = set(Edges)
	Edges = [[Edge[0],Edge[1]] for Edge in Edges]

	# Create an array to store the actual spatial coordinates of all the points
	X = zeros((len(Points),3), float64)
	for i in range(len(Points)):
		X[i,:] = Points[i]

	return X, array(Edges), NewTriangles
	
def MakeTriangle (N, p1, p2, p3, p1p2_EdgeOn = True, p2p3_EdgeOn = True, p3p1_EdgeOn = True, Triangles = None, Points = None):
	"""Makes a triangular sheet and stores it in the array Points (appending values to the end of the array).
	
	p1, p2, p3 specify the corners of the triangle.
	pApB_EdgeOn specifies whether to include the edge between pA and pB."""

	if Points is None:
		Points = []
		TriangleIndexOffset = 0
	else:
		TriangleIndexOffset = len(Points)

	p1, p2, p3 = EnsureNumpyArray(p1), EnsureNumpyArray(p2), EnsureNumpyArray(p3)

	RowOffset = []

	k = 0
	for i in range(N + 1):
		RowOffset.append(k)
		
		t = i / float(N)
		x = (1 - t) * p1 + t * p3
		y = (1 - t) * p2 + t * p3

		if i == 0 and not p1p2_EdgeOn: continue
		for j in range(N + 1 - i):			
			if j == 0 and not p3p1_EdgeOn: continue
			if j == N - i and not p2p3_EdgeOn: continue
			
			if N - i == 0:
				s = 0
			else:
				s = j / float(N - i)			   
			Points.append((1 - s) * x + s * y)

			k += 1

	Points = array(Points)

	if not Triangles is None:
		for i in range(len(RowOffset) - 1):
			Width = RowOffset[i+1] - RowOffset[i]
			for j in range(RowOffset[i], RowOffset[i+1] - 1):
				Triangles.append(Points[[j, j+1, j+Width]])

		for i in range(1, len(RowOffset) - 1):
			Width = RowOffset[i] - RowOffset[i-1]
			for j in range(RowOffset[i], RowOffset[i+1] - 1):
				Triangles.append(Points[[j, j+1, j+1-Width]])

	return Points, Triangles
