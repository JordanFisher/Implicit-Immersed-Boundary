from MatrixUtility import Stack, Unstack
from scipy.sparse.linalg import LinearOperator

from numpy import float64

from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import cg

	
def FluidSolve(fiber, fluid):
	"""Perfom a fluid solve using the actual fluid solver.
	
	The operation includes the spreading of the fiber force and the interpolation back to the fiber."""
	
	fiber.ToGrid(fluid, X = fiber.LaggedX)
	fluid.FluidSolve(fluid.f)
	fiber.FromGrid(fluid, fluid.Output_u, X = fiber.LaggedX)


def Make_KrylovFastEval(fiber, fluid, tree):
	def KrylovFastEval(x):
		fiber.X = x
		fiber.CalcLinearForce()
		tree.FastFluidSolve(fiber)
		
		return x - fluid.dt * fiber.U

	KrylovFastEval.Size = fiber.Size()
		
	return KrylovFastEval
		
def Make_CalcFastResidual(fiber, fluid, tree):
	def CalcFastResidual(x):
		fiber.X = x
		fiber.CalcForceDistribution()
		tree.FastFluidSolve(fiber)
		
		return x - fluid.dt * fiber.U
		
	return CalcFastResidual
	
def Make_KrylovEval(fiber, fluid):
	def KrylovEval(x):
		fiber.X = x
		fiber.CalcLinearForce()
		FluidSolve(fiber, fluid)
		
		return x - fluid.dt * fiber.U

	KrylovEval.Size = fiber.Size()
	
	return KrylovEval

def Make_CalcResidual(fiber, fluid):
	def CalcResidual(x):
		fiber.X = x
		fiber.CalcForceDistribution()
		FluidSolve()
		
		return x - fluid.dt * fiber.U

	return CalcResidual

def StackFunction(f):
	"""Create a stacked version of the function f.
	
	function x -> Stack(f(Unstack(x)))"""
	
	return lambda x: Stack(f(Unstack(x)))

def Operator(f, size=None):
	"""Create a stacked version of the function f, casted as a LinearOperator.
	
	size is the dimension of the operator f."""
	
	if size == None: size = f.Size
	
	operator = LinearOperator((size,size), matvec = StackFunction(f), dtype=float64)
	operator.__call__ = lambda x: operator * x
	
	return operator
	
def StackedCG(f, b, tol = .001, callback=None):
	"""Perform Conjugate Gradient on the system f(x) = b.
	
	The CG iteration is 'Stacked' in the sense that x is converted into a column vector, and f is converted into a matrix.
	The CG code is built in to numpy, but I'm not convinced it's as good as my own CG code."""
	
	f = Operator(f)
	
	dx, info = cg(f, Stack(b), tol=.001, callback=callback)		

	print "CG finished with max residual", abs(f * dx - Stack(b)).max()

	dx = Unstack(dx)

	return dx