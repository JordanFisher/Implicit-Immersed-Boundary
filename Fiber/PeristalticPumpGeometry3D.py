from Tethered import Tethered
import pylab
from numpy import pi, cos, sin, array, zeros, float64, mgrid

from PeristalticPumpGeometry2D import PeristalticPump2D, PumpOrientation, PumpingStyle, PeristalticRadius

class ToSurfaceStyle:
    """Different styles for converting the geometry to a surface."""

    Seamless, Whole, Half = range(3)

class PeristalticPump3D(PeristalticPump2D):
    @staticmethod
    def N_to_dt(N, s):
        """Return the largest stable timestep for an explicit simulation."""

        h = 1. / N
        return .1 * h / s**.5

    def __init__(self, Nb_x, s, Scale = 1., Chi = .1, alpha = pi / 2, Period = 1.,
                 Orientation = PumpOrientation.Horizontal,
                 Style = PumpingStyle.LeftToRight):
        """Set aside the needed memory to handle the ellipse.

        Nb_x is the number of fiber points along a 1D strand in the x direction (Nb = 2 * Nb_x).
        s is a stiffness constant.
        Chi is the channel occlusion ratio. 1 is fully occluded, 0 is fully open.
        Scale is the length of the pump from end to end (the length of one wavelength)"""

        self.Scale = Scale

        self.Chi = float(Chi)
        self.alpha = float(alpha)

        self.Period = float(Period)

        self.Nb_x = Nb_x

        self.Orientation = Orientation

        self.Style = Style

        # Calculate the number of fiber points along a circumference of the perilstaltic pump
        self.Nb_C = int(pi * Nb_x)   

        # Calculate the point separation distance
        self.hb_x, self.hb_C = 1. / self.Nb_x, 2 * pi / self.Nb_C

        Nb = self.Nb_x * self.Nb_C
        Tethered.__init__(self, Nb, Dim = 3, s = s)

        # Initialize x and theta coordinate grids
        self._x, self._theta = zeros(Nb, float64), zeros(Nb, float64)
        for i in range(self.Nb_x):
            x = (i + .5) * self.hb_x

            for j in range(self.Nb_C):
                theta = j * self.hb_C
                
                self._x[i * self.Nb_C + j] = x
                self._theta[i * self.Nb_C + j] = theta
        
        self.SetTime(0, self.X)

    def PeristalticPos(self, x, theta, t):
        """Calculate a point on the pump at time t.

        If Top is True then the point is part of the top curve, otherwise the point is part of the bottom curve."""

        if self.Style == PumpingStyle.LeftToRight:
            t = t / self.Period
        elif self.Style == PumpingStyle.Periodic:
            t = sin(2 * pi * t / self.Period)

        r = PeristalticRadius(x, t, self.Chi, self.alpha)
            
        if hasattr(x, '__iter__'):
            pos = self.X.copy()
            if self.Orientation == PumpOrientation.Horizontal:        
                pos[...,...,0] = x
                pos[...,1] = .5 + r * cos(theta)
                pos[...,2] = .5 + r * sin(theta)
            else:        
                pos[...,0] = .5 + r * sin(theta)                
                pos[...,1] = .5 + r * cos(theta)
                pos[...,2] = x
                
            return pos
        else:
            if self.Orientation == PumpOrientation.Horizontal:        
                return array([x, .5, .5]) + r * array([0, cos(theta), sin(theta)])
            else:        
                return array([.5, .5, x]) + r * array([sin(theta), cos(theta), 0])

    def PeristalticRadius(self, x, t):
        """Calculate the radius of the cross-section through the plane x=x at time t."""

        return PeristalticRadius(x, t, self.Chi, self.alpha)

    def ToSurface(self, Angle = 0, Style = ToSurfaceStyle.Half):
        """Create a 2D array storing the pump's configuration."""

        Nb_x, Nb_C = self.Nb_x, self.Nb_C

        ## Half surface cutaway
        if Style == ToSurfaceStyle.Half:
            offset = int(Nb_C * Angle / 360.0)
            surface = zeros((Nb_x,Nb_C/2,3),float64)

            for i in range(Nb_x):
                surface[i,:] = self.X[i*Nb_C + offset:i*Nb_C + offset + Nb_C/2,:]

        ## Whole surface, with natural cut
        elif Style == ToSurfaceStyle.Whole:
            surface = zeros((Nb_x,Nb_C,3),float64)

            for i in range(Nb_x):
                surface[i,:] = self.X[i*Nb_C:(i+1)*Nb_C,:]

        ## Seamless surface
        elif Style == ToSurfaceStyle.Seamless:
            surface = zeros((Nb_C+1,Nb_x,3),float64)

            for i in range(Nb_C):
                surface[i,:] = self.X[i*Nb_x:(i+1)*Nb_x,:]

            surface[Nb_C,:] = self.X[0:Nb_x,:]

        return surface

    def SetTime(self, t, X=None):
        """Update the configuration of the pump to time t."""

        self.LastSetTime = t
        Tethered.SetTime(self, t)

        if X == None: X = self.Tether

        X[...] = self.PeristalticPos(self._x, self._theta, t)

        #for i in range(self.Nb_x):
            #x = (i + .5) * self.hb_x

            #for j in range(self.Nb_C):
                #theta = j * self.hb_C
                #X[i * self.Nb_C + j] = self.PeristalticPos(x, theta, t)


    def CalcFlux(self, x, fluid):
        """Calculate the fluid flux in the pump past a vertical plane."""

        mask = self.MakeFluidMask_Orthogonal(x, fluid)

        i = int(x / fluid.h[0]) # Get the Eulerian index associated with x

        # Calculate \int_{x^2+y^2<r^2} u_0(x,y,z) dy,
        # where r is the radius of the peristaltic pump at x
        return (fluid.u[i, :, :, 0] * (1-mask)).sum() * fluid.h[1] * fluid.h[2]

    def MakeFluidColorMask_Orthogonal(self, x, fluid, padding=0):
        """Creates a 2D color field that is transparent in the region enclosed by the peristaltic pump.
        The slice is orthogonal to the direction of the pump, the plane (x,?,?).

        padding is an integer that gives the additional number of Eulerian cells above and below the pump to include."""

        N, h = fluid.N, fluid.h

        mask = zeros((N[1],N[2],4))
        t = self.LastSetTime

        # Get the pump radius of the orthognal slice
        r = self.PeristalticRadius(x, t)

        # Calculate the radius from the center [x,.5,.5] for each point in the slice
        grid = mgrid[0:N[1]:1, 0:N[2]:1]
        x, y = grid[0] * h[1], grid[1] * h[2]
        R2 = x**2 + y**2

        # Set the mask to 1 for points with radius greater than the radius of this orthogonal slice
        mask[R2 < r**2] = 1.

        return mask

    def GetVerticalRange(self, fluid, x, padding = 0, t = None):
        """Get two vertical indices that bound the interior of the pump
        on the Eulerian grid along the line (x,?,.5)."""

        if t == None:
            t = self.LastSetTime

        y1 = self.PeristalticPos(x, 0., t)[1]
        y2 = self.PeristalticPos(x, pi, t)[1]

        y1 = int(y1 / fluid.h[1])
        y2 = int(y2 / fluid.h[1])

        y1 = max(0, y1 - padding)
        y2 = min(fluid.N[1], y2 + padding)		

        return y2, y1

    def MakeFluidMask_Orthogonal(self, fluid, padding=0):
        """Creates a 2D scalar field that is 0 in the region enclosed by the peristaltic pump.
        The slice is orthogonal to the direction of the pump, the plane (x,?,?)."""	

        return self.MakeFluidColorMask_Orthogonal(fluid, padding)[...,3]

    def GetSlice(self, theta = 0, X=None):
        """Get a length-wise slice of the pump."""

        if X == None: X = self.X

        i = int(theta / self.hb_C)
        return X[i::self.Nb_C,:]

    def Plot2D(self, X=None, DrawTethers = False, **Scatter_args):
        """Plot a 2D slice slice of the pump."""

        if DrawTethers:
            self.PlotSlice(self.Tether, c='k')
            self.PlotSlice(X, Scatter_args)
            return

        if X == None: X = self.X

        Z = self.GetSlice(0, X=X)
        pylab.scatter(Z[:,0], Z[:,1], **Scatter_args)

        Z = self.GetSlice(pi, X=X)
        pylab.scatter(Z[:,0], Z[:,1], **Scatter_args)

        pylab.xlim(0, 1)
        pylab.ylim(0, 1)