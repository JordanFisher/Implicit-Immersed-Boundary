from FieldUtility import DimOfVectorFieldDomain

def Difference_SecondPeriodic1D(N, h, f, df):
	"""Compute central second order difference of a 1D periodic vector-valued function f.
	
	df must have the same dimensions as f and will store the difference"""
	
	df[1:N-1] = (f[2:] - 2 * f[1:N-1] + f[:N-2]) / h**2
	df[0] = (f[N-1] - 2 * f[0] + f[1]) / h**2
	df[N-1] = (f[N-2] - 2 * f[N-1] + f[0]) / h**2

def Difference_CentralPeriodic(N, h, f, df, axis):
	"""Compute the central difference of a scalar field f,
	in the direction specified by axis, with periodic boundary.
	
	The difference of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2]
	axis must be 'x', 'y', 'z' or 0, 1, 2"""

	if axis == 'x': axis = 0
	if axis == 'y': axis = 1
	if axis == 'z': axis = 2
	
	if axis == 0:
		df[1:N[0]-1,...] = (f[2:,...] - f[:N[0]-2,...]) / (2. * h[0])
		df[N[0]-1,...] = (f[0,...] - f[N[0]-2,...]) / (2. * h[0])
		df[0,...] = (f[1,...] - f[N[0]-1,...]) / (2. * h[0])
	elif axis == 1:
		df[:,1:N[1]-1,...] = (f[:,2:,...] - f[:,:N[1]-2,...]) / (2. * h[1])
		df[:,N[1]-1,...] = (f[:,0,...] - f[:,N[1]-2,...]) / (2. * h[1])
		df[:,0,...] = (f[:,1,...] - f[:,N[1]-1,...]) / (2. * h[1])
	elif axis == 2:
		df[...,1:N[2]-1] = (f[...,2:] - f[...,:N[2]-2]) / (2. * h[2])
		df[...,N[2]-1] = (f[...,0] - f[...,N[2]-2]) / (2. * h[2])
		df[...,0] = (f[...,1] - f[...,N[2]-1]) / (2. * h[2])
	else:
		raise "incorrect axis specification"

def Difference_RightPeriodic(N, h, f, df, axis):
	"""Compute the right difference of a scalar field f,
	df[i] = (df[i+1] - df[i]) / h,
	in the direction specified by axis, with periodic boundary.
	
	The difference of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2]
	axis must be 'x', 'y', 'z' or 0, 1, 2"""

	if axis == 'x': axis = 0
	if axis == 'y': axis = 1
	if axis == 'z': axis = 2
	
	if axis == 0:
		df[:N[0]-1,...] = (f[1:,...] - f[:N[0]-1,...]) / h[0]
		df[N[0]-1,...] = (f[0,...] - f[N[0]-1,...]) / h[0]
	elif axis == 1:
		df[:,:N[1]-1,...] = (f[:,1:,...] - f[:,:N[1]-1,...]) / h[1]
		df[:,N[1]-1,...] = (f[:,0,...] - f[:,N[1]-1,...]) / h[1]
	elif axis == 2:
		df[...,:N[2]-1] = (f[...,1:] - f[...,:N[2]-1]) / h[2]
		df[...,N[2]-1] = (f[...,0] - f[...,N[2]-1]) / h[2]
	else:
		raise "incorrect axis specification"
	
def Difference_LeftPeriodic(N, h, f, df, axis):
	"""Compute the left difference of a scalar field f,
	df[i] = (df[i] - df[i-1]) / h,
	in the direction specified by axis, with periodic boundary.
	
	The difference of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2]
	axis must be 'x', 'y', 'z' or 0, 1, 2"""

	if axis == 'x': axis = 0
	if axis == 'y': axis = 1
	if axis == 'z': axis = 2
	
	if axis == 0:
		df[1:,...] = (f[1:,...] - f[:N[0]-1,...]) / h[0]
		df[0,...] = (f[0,...] - f[N[0]-1,...]) / h[0]
	elif axis == 1:
		df[:,1:,...] = (f[:,1:,...] - f[:,:N[1]-1,...]) / h[1]
		df[:,0,...] = (f[:,0,...] - f[:,N[1]-1,...]) / h[1]
	elif axis == 2:
		df[...,1:] = (f[...,1:] - f[...,:N[2]-1]) / h[2]
		df[...,0] = (f[...,0] - f[...,N[2]-1]) / h[2]
	else:
		raise "incorrect axis specification"
	
def Difference2nd_CentralPeriodic(N, h, f, df, axis):
	"""Compute the 2nd central difference of a scalar field f,
	in the direction specified by axis, with periodic boundary.
	
	The 2nd difference of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2]
	axis must be 'x', 'y', 'z' or 0, 1, 2"""

	if axis == 'x': axis = 0
	if axis == 'y': axis = 1
	if axis == 'z': axis = 2
	
	h2 = h**2
	
	if axis == 0:
		df[1:N[0]-1,...] = f[2:,...] + f[:N[0]-2,...]
		df[N[0]-1,...] =   f[0,...]  + f[N[0]-2,...]
		df[0,...] =        f[1,...]  + f[N[0]-1,...]
		df[...] = (df - 2. * f) / h2[0]
	elif axis == 1:
		df[:,1:N[1]-1,...] = f[:,2:,...] + f[:,:N[1]-2,...]
		df[:,N[1]-1,...]   = f[:,0,...]  + f[:,N[1]-2,...]
		df[:,0,...]        = f[:,1,...]  + f[:,N[1]-1,...]
		df[...] = (df - 2. * f) / h2[1]
	elif axis == 2:
		df[...,1:N[2]-1] = f[...,2:] + f[...,:N[2]-2]
		df[...,N[2]-1]   = f[...,0]  + f[...,N[2]-2]
		df[...,0]        = f[...,1]  + f[...,N[2]-1]
		df[...] = (df - 2. * f) / h2[2]
	else:
		raise "incorrect axis specification"

def Laplacian(N, h, f, df):
	"""Compute the Laplacian of a scalar field f, with periodic boundary.
	
	The Laplacian of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2]"""

	Dim = f.ndim
	df *= 0
	
	for axis in range(Dim):
		derivative = df * 0
		Difference2nd_CentralPeriodic(N, h, f, derivative, axis)
		df += derivative

def VectorLaplacian(N, h, f, df=None):
	"""Compute the Laplacian of a vector field f, with periodic boundary.
	
	The Laplacian of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2] x n, for some n"""

	if df == None:
		df = 0 * f
	
	Dim = DimOfVectorFieldDomain(f)
	df *= 0
	
	n = f.shape[Dim] # The number of components each location in the field has
	
	for i in range(n):
		Laplacian(N, h, f[...,i], df[...,i])
	
	return df

class DifferenceTypes:
	CentralPeriodic, LeftPeriodic, RightPeriodic  = range(3)

def Difference(N, h, f, df, axis, DifferenceType = DifferenceTypes.CentralPeriodic):
	"""Compute a difference of a scalar field f,
	in the direction specified by axis, with periodic boundary.
	
	The difference of f will be stored in df
	N and h must be arrays of length 2 or 3
	f, df must be arrays with dimension N[0]xN[1]xN[2]
	axis must be 'x', 'y', 'z' or 0, 1, 2"""
	
	if DifferenceType == DifferenceTypes.CentralPeriodic:
		return Difference_CentralPeriodic(N, h, f, df, axis)
	elif DifferenceType == DifferenceTypes.LeftPeriodic:
		return Difference_LeftPeriodic(N, h, f, df, axis)
	elif DifferenceType == DifferenceTypes.RightPeriodic:
		return Difference_RightPeriodic(N, h, f, df, axis)
	else:
		raise "Unknown difference scheme requested."
	
def TotalDifference(N, h, u, du, DifferenceType = DifferenceTypes.CentralPeriodic):
	"""Take a 2D vector field or 3D vector field u and store a 2x2 matrix field or a 3x3 matrix field of differences, df.
	
	N and h must be arrays of length 2 or 3, refered to as dim.
	f must be an array with dimensions N[0]xN[1]xN[2]x n, for any positive n
	df must be an array with dimension N[0]xN[1]xN[2]x n x dim"""

	Dim = DimOfVectorFieldDomain(u)

	n = u.shape[Dim] # The number of components each location in the field has
	for axis in range(Dim): # Each direction to derive in: 'x', 'y', 'z'
		for i in range(n):
			Difference(N, h, u[...,i], du[...,i,axis], axis, DifferenceType)