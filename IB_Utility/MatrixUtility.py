from numpy import zeros, dtype, float64, eye

def DssMatrix3D (Nb):
    """Computes the matrix of the operator Dss, the periodic second difference."""

    Nb2 = Nb**2

    D = zeros((3*Nb,3*Nb), dtype=float64)
    D[:Nb,:Nb] -= 2 * eye(Nb,Nb) * Nb2
    D[1:Nb,:Nb-1] += eye(Nb-1,Nb-1) * Nb2
    D[:Nb-1,1:Nb] += eye(Nb-1,Nb-1) * Nb2
    D[Nb-1,0] = D[0,Nb-1] = Nb2

    D[Nb:2*Nb,Nb:2*Nb] = D[:Nb,:Nb]
    D[2*Nb:,2*Nb:] = D[:Nb,:Nb]

    return D

def Stack(A, _A = 0):
    """Takes a vector array A and returns a scalar array.
    If A = (X,Y,Z), where X,Y,Z are scalar arrays then
    _A is the concatenation of the three arrays."""
    
    if len(A.shape) == 2:
        N, d = A.shape
        if type(_A) == type(0):
            _A = zeros(d * N, A.dtype)
        for i in range(d):
            _A[i*N:(i+1)*N] = A[:,i]

    elif len(A.shape) == 4:
        N, M, d1, d2 = A.shape
        if N != M or d1 != d2:
            raise "A is not square or has incorrect dimensions"
        
        if type(_A) == type(0):
            _A = zeros((d1*N,d2*M),float64)
            
        for i in range(d1):
            for j in range(d2):
                _A[i*N:(i+1)*N,j*N:(j+1)*N] = A[:,:,i,j]

    return _A        

def Unstack(A, Dim = 3):
    """Takes a scalar array A and returns the 3-vector array _A that it represents.
	If Dim = 2 then a 2-vectory array is returned."""

    l = len(A)
    N = l / Dim
    _A = zeros((N, Dim), A.dtype)
    for i in range(Dim):
        _A[:,i] = A[i*N:(i+1)*N]

    return _A
    
def ComputeExactMatrix (fiber, fluid):
    """For a given fiber/fluid configuration computes the matrix represenation
    in stacked form of the fiber force to fiber motion operator M_n"""
    
    A = zeros((fiber.Nb, fiber.Nb,3,3), float64)

    holdu = fluid.u.copy()
    holdF = fiber.F.copy()
    holdU = fiber.U.copy()

    for j in range(fiber.Nb):
        for l in range(3):
            fiber.F *= 0
            fiber.F[j][l] = 1.
            fluid.u *= 0

            FiberToGrid (fluid.N, fluid.h, fiber.Nb, fiber.hb, fiber.X, fiber.F, fluid.f)
            fluid.FluidSolve(fluid.f)
            GridToFiber (fluid.N, fluid.h, fiber.Nb, fiber.hb, fiber.X, fluid.u, fiber.U)

            A[:,j,l,:] = fiber.U.copy()

    fluid.u = holdu
    fiber.F = holdF
    fiber.U = holdU

    return A
