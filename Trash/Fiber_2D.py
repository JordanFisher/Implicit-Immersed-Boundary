from numpy import zeros
from IBUtility import EnsureNumpyArray
import Fiber, Tethered

def TriangulateTriangle(Triangle):
    Triangles = []
    a, b, c = Triangle
    ab, bc, ac = (a+b)/2, (b+c)/2, (a+c)/2
    return [(a, ab, ac), (b, ab, bc), (c, bc, ac), (ab, bc, ac)]

def Triangulate(Triangles):
    return [NewTriangle for Triangle in Triangles for NewTriangle in TriangulateTriangle(Triangle)]

##def TrianglesToIndices(Triangles, X):
##    IndexedTriangles = []
##    for Triangle in Triangles:
##        IndexedTriangles.append([X.Triangle[0]

def TrianglesToPoints(Triangles):
    Points = []
    for Triangle in Triangles:
        Points.append((Triangle[0][0], Triangle[0][1], Triangle[0][2]))
        Points.append((Triangle[1][0], Triangle[1][1], Triangle[1][2]))
        Points.append((Triangle[2][0], Triangle[2][1], Triangle[2][2]))        

    PointDic = {}
    map(lambda x: PointDic.__setitem__(x, 0), Points)
    Points = PointDic.keys()
    i = 0
    for key in Points:
        PointDic[key] = i
        i += 1

    NewTriangles = [[PointDic[(Triangle[i][0], Triangle[i][1], Triangle[i][2])] for i in range(3)] for Triangle in Triangles]

    Edges = []
    for Triangle in NewTriangles:
        for i in range(-1,2):
            if Triangle[i] < Triangle[i+1]:
                Edges.append( (Triangle[i], Triangle[i+1]) )
            else:
                Edges.append( (Triangle[i+1], Triangle[i]) )
    Edges = set(Edges)
    Edges = [[Edge[0],Edge[1]] for Edge in Edges]

    X = zeros((len(Points),3), float64)
    for i in range(len(Points)):
        X[i,:] = Points[i]

    return X, NewTriangles, array(Edges)

def ProjectPointToSphere(p, c, r):
    dif = p - c
    l = dot(dif, dif)**.5
    return dif * r / l + c

def ProjectToSphere(Triangles, c, r):
    for Triangle in Triangles:
        Triangle[0][:] = ProjectPointToSphere(Triangle[0], c, r)
        Triangle[1][:] = ProjectPointToSphere(Triangle[1], c, r)
        Triangle[2][:] = ProjectPointToSphere(Triangle[2], c, r)

def MakeSphere(c, r):
    Triangles = []
##    Triangles.append((array([0,0,r],float64)+c, array([r,0,0],float64)+c, array([0,r,0],float64)+c))
##    Triangles.append((array([0,0,r],float64)+c, array([r,0,0],float64)+c, array([0,-r,0],float64)+c))
##    Triangles.append((array([0,0,r],float64)+c, array([-r,0,0],float64)+c, array([0,r,0],float64)+c))
##    Triangles.append((array([0,0,r],float64)+c, array([-r,0,0],float64)+c, array([0,-r,0],float64)+c))

    phi = (1 + 5**.5) / 2
    R = (1 + phi**2)**.5 # Radius of standard icosahedron
    ratio = r / R # Amount to scale the icosahedron so it has radius r
    phi *= ratio
    
    Triangles.append((array([ratio,0,phi],float64)+c, array([-ratio,0,phi],float64)+c, array([0,phi,ratio],float64)+c))
    Triangles.append((array([ratio,0,phi],float64)+c, array([-ratio,0,phi],float64)+c, array([0,-phi,ratio],float64)+c))    

    Triangles.append((array([-phi,ratio,0],float64)+c, array([-ratio,0,phi],float64)+c, array([0,phi,ratio],float64)+c))
    Triangles.append((array([-phi,-ratio,0],float64)+c, array([-ratio,0,phi],float64)+c, array([0,-phi,ratio],float64)+c))
    Triangles.append((array([-phi,-ratio,0],float64)+c, array([-ratio,0,phi],float64)+c, array([-phi,ratio,0],float64)+c))    

    Triangles.append((array([phi,ratio,0],float64)+c, array([ratio,0,phi],float64)+c, array([0,phi,ratio],float64)+c))
    Triangles.append((array([phi,-ratio,0],float64)+c, array([ratio,0,phi],float64)+c, array([0,-phi,ratio],float64)+c))
    Triangles.append((array([phi,-ratio,0],float64)+c, array([ratio,0,phi],float64)+c, array([phi,ratio,0],float64)+c))    

    Triangles.append((array([phi,-ratio,0],float64)+c, array([0,-phi,ratio],float64)+c, array([0,-phi,-ratio],float64)+c))
    Triangles.append((array([-phi,-ratio,0],float64)+c, array([0,-phi,ratio],float64)+c, array([0,-phi,-ratio],float64)+c))

    Triangles.append((array([phi,ratio,0],float64)+c, array([0,phi,ratio],float64)+c, array([0,phi,-ratio],float64)+c))
    Triangles.append((array([-phi,ratio,0],float64)+c, array([0,phi,ratio],float64)+c, array([0,phi,-ratio],float64)+c))


    Triangles.append((array([ratio,0,-phi],float64)+c, array([-ratio,0,-phi],float64)+c, array([0,phi,-ratio],float64)+c))
    Triangles.append((array([ratio,0,-phi],float64)+c, array([-ratio,0,-phi],float64)+c, array([0,-phi,-ratio],float64)+c))    

    Triangles.append((array([-phi,ratio,0],float64)+c, array([-ratio,0,-phi],float64)+c, array([0,phi,-ratio],float64)+c))
    Triangles.append((array([-phi,-ratio,0],float64)+c, array([-ratio,0,-phi],float64)+c, array([0,-phi,-ratio],float64)+c))
    Triangles.append((array([-phi,-ratio,0],float64)+c, array([-ratio,0,-phi],float64)+c, array([-phi,-ratio,0],float64)+c))    

    Triangles.append((array([phi,ratio,0],float64)+c, array([ratio,0,-phi],float64)+c, array([0,phi,-ratio],float64)+c))
    Triangles.append((array([phi,-ratio,0],float64)+c, array([ratio,0,-phi],float64)+c, array([0,-phi,-ratio],float64)+c))
    Triangles.append((array([phi,-ratio,0],float64)+c, array([ratio,0,-phi],float64)+c, array([phi,-ratio,0],float64)+c))    

    for i in range(2):
        Triangles = Triangulate(Triangles)
        ProjectToSphere(Triangles, c, r)

    ProjectToSphere(Triangles, c, r)
    return TrianglesToPoints(Triangles)

def MakeTriangle (N, p1, p2, p3, p1p2_EdgeOn = True, p2p3_EdgeOn = True, p3p1_EdgeOn = True, Triangles = None, Points = None):
    """Makes a triangular sheet and stores it in Points
    p1, p2, p3 specify the corners of the triangle"""

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
                #Triangles.append(array([j, j+1, j+Width]) + TriangleIndexOffset)
                Triangles.append(Points[[j, j+1, j+Width]])

        for i in range(1, len(RowOffset) - 1):
            Width = RowOffset[i] - RowOffset[i-1]
            for j in range(RowOffset[i], RowOffset[i+1] - 1):
                #Triangles.append(array([j, j+1, j+1-Width]) + TriangleIndexOffset)
                Triangles.append(Points[[j, j+1, j+1-Width]])

    return Points, Triangles



def MakeSphere2 (N, c, r):
    c = EnsureNumpyArray(c)

    Triangles, Points = [], []
##    MakeTriangle (N, c + [0,-r,0], c + [r,0,0], c + [0,0,r], Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,-r,0], c + [-r,0,0], c + [0,0,r], p3p1_EdgeOn = False, Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,r,0], c + [r,0,0], c + [0,0,r], p2p3_EdgeOn = False, Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,r,0], c + [-r,0,0], c + [0,0,r], p3p1_EdgeOn = False, p2p3_EdgeOn = False, Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,-r,0], c + [r,0,0], c + [0,0,-r], p1p2_EdgeOn = False, Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,-r,0], c + [-r,0,0], c + [0,0,-r], p3p1_EdgeOn = False, p1p2_EdgeOn = False, Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,r,0], c + [r,0,0], c + [0,0,-r], p2p3_EdgeOn = False, p1p2_EdgeOn = False, Triangles = Triangles, Points = Points)
##    MakeTriangle (N, c + [0,r,0], c + [-r,0,0], c + [0,0,-r], p3p1_EdgeOn = False, p2p3_EdgeOn = False, p1p2_EdgeOn = False, Triangles = Triangles, Points = Points)
    MakeTriangle (N, c + [0,-r,0], c + [r,0,0], c + [0,0,r], Triangles = Triangles)
    MakeTriangle (N, c + [0,-r,0], c + [-r,0,0], c + [0,0,r], Triangles = Triangles)
    MakeTriangle (N, c + [0,r,0], c + [r,0,0], c + [0,0,r], Triangles = Triangles)
    MakeTriangle (N, c + [0,r,0], c + [-r,0,0], c + [0,0,r], Triangles = Triangles)
    MakeTriangle (N, c + [0,-r,0], c + [r,0,0], c + [0,0,-r], Triangles = Triangles)
    MakeTriangle (N, c + [0,-r,0], c + [-r,0,0], c + [0,0,-r], Triangles = Triangles)
    MakeTriangle (N, c + [0,r,0], c + [r,0,0], c + [0,0,-r], Triangles = Triangles)
    MakeTriangle (N, c + [0,r,0], c + [-r,0,0], c + [0,0,-r], Triangles = Triangles)

    Points, Triangles, Edges = TrianglesToPoints(Triangles)

    Points = array(Points)

    dif = Points - c
    difsq = dif**2
    l = (difsq[:,0] + difsq[:,1] + difsq[:,2])**.5
    dif[:,0] *= r / l
    dif[:,1] *= r / l
    dif[:,2] *= r / l

    Points = dif + c

    return Points, array(Triangles), array(Edges)

        
class Fiber_2D(Tethered):
    """A 1D fiber"""
    
    def __init__ (self, Nb = 0, s = 0, Points = None, Links = None, Tether = None):
        """Nb gives the number of points used to describe the fiber
        s is an elasticity constant"""

        if not Points is None:
            Nb = Points.shape[0]
            self.X = Points
        else:        
            self.X = zeros((Nb,3),float64)
        self.Tether = Tether
        self.Links = Links
        self.F = zeros((Nb,3),float64)
        self.U = zeros((Nb,3),float64)
        self.Nb = Nb
        self.N = int(Nb**.5)
        self.hb2 = self.hb = 1. / self.Nb
        self.s = s

    def CalcForceDistribution (self):
        """Calculate the force distribution along the fiber"""

        self.F *= 0

        if not self.Tether is None:
            self.F += self.s * (self.Tether - self.X)
            
        if not self.Links is None:
            if self.Links_s is None:
                force = self.s * (self.X[self.Links[:,1]] - self.X[self.Links[:,0]])
            else:
                force = self.X[self.Links[:,1]] - self.X[self.Links[:,0]]
                force[:,0] *= self.Links_s
                force[:,1] *= self.Links_s
                force[:,2] *= self.Links_s                
            for i in range(3):
                binc = bincount(self.Links[:,0], force[:,i])
                self.F[:binc.shape[0], i] += binc
                binc = bincount(self.Links[:,1], force[:,i])
                self.F[:binc.shape[0], i] -= binc

##                import code; code.interact(local=dict(globals(), **locals()))

    def MakeSquare (self, N, p1, p2, p3, h):
        """Make the fiber into a square sheet
        p1, p2, p3 are sequential corners of the square"""
        
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

    def SaveState(self, Saver):
        Saver.Save(self.X, 'X_' + str(Saver.CurStep))
