from numpy import array, zeros, float64, complex64, fft, real, dot

def Dot(a, b):
	"""Calculate the dot product of a * b.
	
	a, b are arrays of vectors"""
	
	return (a * b).sum()

def FieldDot(a, b, c=None, Dim=3):
	"""Calculate a field c defined via c[j,k,l] = a[j,k,l] * b[j,k,l].
	
	a, b may be scalar fields, vector fields or Dim x Dim matrix fields, where Dim is the number of spatial dimensions of the domain."""
	
	if Dim == 2:
		return FieldDot2D(a, b, c)
	elif Dim == 3:
		return FieldDot3D(a, b, c)
	else:
		raise Exception("field is neither a 2D field nor a 3D field")

def FieldDot2D(a, b, c=None):
	"""Calculate a field c defined via c[j,k,l] = a[j,k,l] * b[j,k,l].
	
	a, b may be scalar fields, vector fields or 3x3 matrix fields."""
	
	ScalarDim, VectorDim, MatrixDim = 2, 3, 4
	
	if a.ndim == MatrixDim and b.ndim == VectorDim:
		# matrix vector multiplication
		if c == None: c = zeros(b.shape, a.dtype)
		
		c[...,0] = a[...,0,0] * b[...,0] + a[...,0,1] * b[...,1]
		c[...,1] = a[...,1,0] * b[...,0] + a[...,1,1] * b[...,1]
	elif a.ndim == VectorDim and b.ndim == VectorDim:
		# vector vector dot product
		if c == None: c = zeros((a.shape[0], a.shape[1]), a.dtype)
		
		c[...] = a[...,0] * b[...,0] + a[...,1] * b[...,1]
	elif a.ndim == VectorDim and b.ndim == ScalarDim:
		# vector scalar multiplication
		if c == None: c = zeros(a.shape, a.dtype)
		
		c[...,0] = a[...,0] * b[...]
		c[...,1] = a[...,1] * b[...]
	elif a.ndim == ScalarDim and b.ndim == VectorDim:
		# scalar vector multiplication
		if c == None: c = zeros(b.shape, a.dtype)
		
		c[...,0] = a[...] * b[...,0]
		c[...,1] = a[...] * b[...,1]
	else:
		raise Exception("DotField error: code has not yet defined this operation")
	
	return c
		
def FieldDot3D(a, b, c=None):
	"""Calculate a field c defined via c[j,k,l] = a[j,k,l] * b[j,k,l].
	
	a, b may be scalar fields, vector fields or 3x3 matrix fields."""
	
	ScalarDim, VectorDim, MatrixDim = 3, 4, 5
	
	if a.ndim == MatrixDim and b.ndim == VectorDim:
		# matrix vector multiplication
		if c == None: c = zeros(b.shape, a.dtype)
		
		c[...,0] = a[...,0,0] * b[...,0] + a[...,0,1] * b[...,1] + a[...,0,2] * b[...,2]
		c[...,1] = a[...,1,0] * b[...,0] + a[...,1,1] * b[...,1] + a[...,1,2] * b[...,2]
		c[...,2] = a[...,2,0] * b[...,0] + a[...,2,1] * b[...,1] + a[...,2,2] * b[...,2]
	elif a.ndim == VectorDim and b.ndim == VectorDim:
		# vector vector dot product
		if c == None: c = zeros((a.shape[0], a.shape[1], a.shape[2]), a.dtype)
		
		c[...] = a[...,0] * b[...,0] + a[...,1] * b[...,1] + a[...,2] * b[...,2]
	elif a.ndim == VectorDim and b.ndim == ScalarDim:
		# vector scalar multiplication
		if c == None: c = zeros(a.shape, a.dtype)
		
		c[...,0] = a[...,0] * b[...]
		c[...,1] = a[...,1] * b[...]
		c[...,2] = a[...,2] * b[...]
	elif a.ndim == ScalarDim and b.ndim == VectorDim:
		# scalar vector multiplication
		if c == None: c = zeros(b.shape, a.dtype)
		
		c[...,0] = a[...] * b[...,0]
		c[...,1] = a[...] * b[...,1]
		c[...,2] = a[...] * b[...,2]
	else:
		raise Exception("DotField error: code has not yet defined this operation")
	
	return c

def ComponentSupNorm(a):
	"""Calculate the sup norm of each component of the field a."""
	
	Dim = FieldVectorDim(a)
	
	sup = zeros(Dim)
	for i in range(Dim):
		sup[i] = abs(a[...,i]).max()
	
	return sup
	
def L2(a):
	"""Calculate the l2 norm of a.
	
	a must be a vector field."""
	
	return Dot(a, a)**.5
	
def FieldVectorDim(a):
	"""Dimension of the vectors comprising a field a. For example, if a is a 100x100x3 vector field, 3 is returned."""

	shape = a.shape
	Dim = a.shape[len(shape) - 1]
	
	return Dim
	
def FieldAverage(a):
	"""Calculate the average value of a vector field a."""
	
	shape = a.shape
	Dim = FieldVectorDim(a)

	Elements = 1
	for i in range(len(shape) - 1):
		Elements *= shape[i]

	average = zeros(Dim)
	
	for i in range(Dim):
		average[i] = a[...,i].sum() / Elements
		
	return average

def DimOfVectorFieldDomain(f):
	"""Return the dimension of the domain of the vector field f."""
	
	if f.ndim == 3:
		return 2
	elif f.ndim == 4:
		return 3
	else:
		raise Exception("field is neither a 2D vector field nor a 3D vector field")
	
def FieldFFT(a, a_hat):
	"""Calculate the component-wise 2D or 3D FFT of a vector field a, and stores it in a_hat."""
	
	if DimOfVectorFieldDomain(a) == 2:
		a_hat[:,:,0], a_hat[:,:,1] = fft.fftn(a[:,:,0]), fft.fftn(a[:,:,1])
	else:
		a_hat[...,0], a_hat[...,1], a_hat[...,2] = fft.fftn(a[...,0]), fft.fftn(a[...,1]), fft.fftn(a[...,2])

	return a_hat


def FieldIFFT(a_hat, a):
	"""Calculate the component-wise 2D or 3D inverse FFT of a vector field a_hat, and stores it in a."""

	if DimOfVectorFieldDomain(a) == 2:
		a[:,:,0], a[:,:,1] = real(fft.ifftn(a_hat[:,:,0])), real(fft.ifftn(a_hat[:,:,1]))
	else:
		a[...,0], a[...,1], a[...,2] = real(fft.ifftn(a_hat[...,0])), real(fft.ifftn(a_hat[...,1])), real(fft.ifftn(a_hat[...,2]))

	return a
