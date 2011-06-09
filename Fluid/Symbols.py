from numpy import mgrid, sin, cos, pi, float64, complex64, zeros

def InitSymbols (N, h):
	"""Initialize discrete Fourier space symbols.
	
	Returns grid, WLambda_Hat, SLambda_Hat, D_Hat.

	N and h must be arrays of length 2 or 3."""

	Dim = len(N) # The number of spatial dimensions
	
	if Dim == 2:
		grid = mgrid[0:N[0]:1, 0:N[1]:1]
	elif Dim == 3:
		grid = mgrid[0:N[0]:1, 0:N[1]:1, 0:N[2]:1]
	else:
		raise Exception("Domain must be either 2D or 3D")
	
	WLambda_Hat = InitWideLaplacian (N, h, grid)
	SLambda_Hat = InitShortLaplacian (N, h, grid)
	D_Hat = InitDSymbol (N, h, grid)

	# Create a class to hold the various symbols
	class SymbolContainer:
		pass
		
	a = SymbolContainer()
	
	a.WLambda_Hat, a.SLambda_Hat = WLambda_Hat, SLambda_Hat
	a.D_Hat = D_Hat

	return a


def InitWideLaplacian (N, h, grid):
	""""Initialize the wide Laplacian Fourier symbol via the formula.
	
	In 3D the formula is
	
	LambdaHat[j,k,l] = -sin(2*pi*j / N[0])**2 / h[0]**2 - sin(2*pi*k / N[1])**2 / h[1]**2 - sin(2*pi*l / N[2])**2 / h[2]**2
	
	N and h must be arrays of length 2 or 3
	grid must be integer arrays of dimension N[0]xN[1] or N[0]xN[1]xN[2]"""
	
	Dim = len(N) # The number of spatial dimensions
	
	if Dim == 2:
		Lambda = -sin(2*pi*grid[0] / N[0])**2 / h[0]**2 - sin(2*pi*grid[1] / N[1])**2 / h[1]**2
		Lambda[0,0] = 1.
	else:
		Lambda = -sin(2*pi*grid[0] / N[0])**2 / h[0]**2 - sin(2*pi*grid[1] / N[1])**2 / h[1]**2 - sin(2*pi*grid[2] / N[2])**2 / h[2]**2
		Lambda[0,0,0] = 1.
	
	return Lambda

def InitShortLaplacian (N, h, grid):
	""""Initialize the short Laplacian Fourier symbol via the formula.
	
	In 3D the formula is
	
	LambdaHat[j,k,l] = (2 * cos(2.*pi*j / N[0]) - 2.) / h[0]**2 + (2 * cos(2.*pi*k / N[1]) - 2.) / h[1]**2 + (2 * cos(2.*pi*l / N[2]) - 2.) / h[2]**2
	
	N and h must be arrays of length 2 or 3
	grid must be integer arrays of dimension N[0]xN[1] or N[0]xN[1]xN[2]"""
	
	Dim = len(N) # The number of spatial dimensions
	
	if Dim == 2:
		Lambda = (2 * cos(2.*pi*grid[0] / N[0]) - 2.) / h[0]**2 + (2 * cos(2.*pi*grid[1] / N[1]) - 2.) / h[1]**2
		Lambda[0,0] = 1.
	else:
		Lambda = (2 * cos(2.*pi*grid[0] / N[0]) - 2.) / h[0]**2 + (2 * cos(2.*pi*grid[1] / N[1]) - 2.) / h[1]**2 + (2 * cos(2.*pi*grid[2] / N[2]) - 2.) / h[2]**2
		Lambda[0,0,0] = 1.

	return Lambda


def InitDSymbol (N, h, grid):
	""""Initialize the Fourier symbol for a total difference operator.
	
	Returns a vector field of dimensions N[0]xN[1] or N[0]xN[1]xN[2]
	
	N and h must be arrays of length 3
	grid must be integer arrays of dimension N[0]xN[1] or N[0]xN[1]xN[2]"""

	Dim = len(N) # The number of spatial dimensions

	if Dim == 2:
		D_Hat = zeros((N[0],N[1],2),complex64)
	else:
		D_Hat = zeros((N[0],N[1],N[2],3),complex64)
	
	for i in range(Dim):
		D_Hat[...,i] = InitCentralDSymbol(N, h, grid, i)

	return D_Hat

def InitCentralDSymbol (N, h, grid, axis):
	""""Initialize the Fourier symbol for a central derivative in x or y or z.
	For x the formula is
	
	DxHat[j,k,l] = 1j * sin(2*pi*j / N[0]) / h[0]
	
	N and h must be arrays of length 3
	grid must be integer arrays of dimension N[0]xN[1] or N[0]xN[1]xN[2]"""

	return 1j * sin(2 * pi * grid[axis] / N[axis]) / h[axis]
