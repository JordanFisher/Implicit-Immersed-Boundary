from Tethered import Tethered
import pylab as m
from numpy import pi, cos, sin, array, zeros, float64, int32

def JaffrinMeanFlow(chi, alpha):
    """Calculate the mean flow for a channel with aspect ratio alpha and occlusion ratio chi
    via Jaffrin's assymptotic formula."""

    a = alpha
    x = chi

    if hasattr(x, '__iter__'):
        x = m.array(x)    
    
    num = 15. * x**2 + 2. * a**2 * (4. * (1. - x**2)**2.5 + (7. * x**2 - 4.) * (1. - x**2))
    denom = x * (5. * (2. + x**2) + 6. * a**2 * x**2 * (1 - x**2))

    return num / denom

def CalculateMeanFlow(flow, chi, alpha):
    """Calculate the normalized mean flow using the formula
    \pi / (\alpha \chi) \Phi,
    where \Phi is the non-normalized mean flow."""

    if hasattr(flow, '__iter__'):
        return list(-m.pi * m.array(flow) / (alpha * m.array(chi_values)))
    else:
        return -m.pi * flow / (alpha * chi)

class PumpingStyle:
    """Different choices for how the pump moves."""

    LeftToRight, Periodic = range(2)

class PumpOrientation:
    """Orientation for a peristaltic pump."""

    Horizontal, Vertical = range(2)

def PeristalticRadius(x, t, Chi, alpha = .25):
    """Calculate the radius of the cross-section through the plane x=x at time t."""

    return (alpha / (2. * pi)) * (1 + Chi * sin(2*pi*(x - t)))

def PeristalticPos(x, t, Chi, alpha, Top, Orientation = PumpOrientation.Horizontal):
    """Calculate a point on the pump at time t.

    If Top is True then the point is part of the top curve, otherwise the point is part of the bottom curve."""

    r = PeristalticRadius(x, t, Chi, alpha)

    if Orientation == PumpOrientation.Horizontal:
        return array([x, .5 + r * (1 if Top else -1)])
    else:
        return array([.5 + r * (1 if Top else -1), x])

class PeristalticPump2D(Tethered):
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
        self.alpha = alpha

        self.Period = float(Period)

        self.Nb_x = Nb_x

        self.Orientation = Orientation

        self.Style = Style

        # Calculate the point separation distance
        self.hb_x = 1. / self.Nb_x

        Tethered.__init__(self, 2 * self.Nb_x, Dim = 2, s = s)

        self.SetTime(0, self.X)

    def PeristalticPos(self, x, t, Top):
        """Calculate a point on the pump at time t.

        If Top is True then the point is part of the top curve, otherwise the point is part of the bottom curve."""

        if self.Style == PumpingStyle.LeftToRight:
            t = t / self.Period
        elif self.Style == PumpingStyle.Periodic:
            t = sin(2 * pi * t / self.Period)

        return self.Scale * PeristalticPos(x / self.Scale, t, self.Chi, self.alpha, Top, self.Orientation)

    def CalcFlux(self, x, fluid):
        """Calculate the fluid flux in the pump past a vertical line."""

        y1, y2 = self.GetVerticalRange(fluid, x, 0)

        i = int(x / fluid.h[0]) # Get the Eulerian index associated with x

        # Calculate \int_y1^y2 u_0(x,y) dy
        return fluid.u[i, y1:y2, 0].sum() * fluid.h[1]

    def StartTrackingFlowRate(self, fluid, steps):
        """Start tracking the flow rate of the pump.
        steps is an integer and reprents how many timesteps to wait between updating."""

        self.steps = steps

        # Used to store the flux at every x-value
        #self.CurrentFlux = zeros(fluid.N[0], float64)

        self.Flux = []
        self.MeanFlow = []
        self.NormalizedFlux = []
        self.NormalizedMeanFlow = []

    def UpdateFlowRate(self, fluid, Time):
        """Update the mean flow calculation.
        Store the new mean flow in self.MeanFlow."""

        if Time.CurStep % self.steps != 0:
            return

        # Recalculate the flux at every x-value
        #for i in range(fluid.N[0]):
            #x = i * fluid.h[0]			
            #self.CurrentFlux[i] = self.CalcFlux(x, fluid)

        # Calculate the flux at the middle of the domain
        x = fluid.dims[0] / 2.
        flux = self.CalcFlux(x, fluid)
        self.Flux.append(flux)

        # Calculate the normalized flux
        normalized_flux = CalculateMeanFlow(flux, self.Chi, self.alpha)		
        self.NormalizedFlux.append(normalized_flux)

        Period = abs(self.Period)
        if Time.t > abs(Period):
            # How many frames to include in the mean flow calculation
            # Should include all frames from the most recent period of the pump
            frames = int(abs(Period) / (Time.dt * self.steps))

            # Average over the last period of the pump
            mean = sum(self.Flux[-frames:]) / frames
            normalized_mean = sum(self.NormalizedFlux[-frames:]) / frames

            self.MeanFlow.append(mean)
            self.NormalizedMeanFlow.append(normalized_mean)

    def MakeFluidColorMask(self, fluid, padding=0):
        """Creates a color field that is transparent in the region enclosed by the peristaltic pump.

        padding is an integer that gives the additional number of Eulerian cells above and below the pump to include."""

        mask = zeros((fluid.N[0],fluid.N[1],4))

        for i in range(fluid.N[0]):
            x = i * fluid.h[0]

            y1, y2 = self.GetVerticalRange(fluid, x, padding)

            if self.Orientation == PumpOrientation.Horizontal:
                mask[i,:y1+1,:] += 1
                mask[i,y2:,:] += 1 
            else:
                mask[:y1+1,i,:] += 1
                mask[y2:,i,:] += 1 

        return mask

    def GetVerticalRange(self, fluid, x, padding = 0, t = None):
        """Get two vertical indices that bound the interior of the pump
        on the Eulerian grid along the vertical line at x."""

        if t == None:
            t = self.LastSetTime

        y1 = self.PeristalticPos(x, t, Top = False)[1]
        y2 = self.PeristalticPos(x, t, Top = True)[1]

        y1 = int(y1 / fluid.h[1])
        y2 = int(y2 / fluid.h[1])

        y1 = max(0, y1 - padding)
        y2 = min(fluid.N[1], y2 + padding)		

        return y1, y2

    def MakeFluidMask(self, fluid, padding=0):
        """Creates a scalar field that is 0 in the region enclosed by the peristaltic pump."""	

        return self.MakeFluidColorMask(fluid, padding)[...,3]

    def SetTime(self, t, X=None):
        """Update the configuration of the pump to time t."""

        self.LastSetTime = t
        Tethered.SetTime(self, t)

        if X == None: X = self.Tether

        for i in range(self.Nb_x):
            x = self.Scale * (i + .5) * self.hb_x

            X[i] = self.PeristalticPos(x, t, Top = True)
            X[i + self.Nb_x] = self.PeristalticPos(x, t, Top = False)

    def MakeModFuncs(self, fluid, Mask):
        if Mask:
            mask = self.MakeFluidColorMask(fluid).swapaxes(0,1)	

        def AdditionalPlot():
            if Mask:
                m.imshow(mask, interpolation='bilinear', extent=fluid.PlotExtent2D())
            self.Plot2D(s = 1)
            fluid.SetExtent()

        def ModVal(val):
            shape = val.shape
            middleVal = val[shape[0]/2,shape[1]/2]

            if Mask:
                return middleVal * mask[...,3].T + val * (1 - mask[...,3].T)
            else:
                return val

        return AdditionalPlot, ModVal

    def VectorPlot(self, fluid, Mask=True):
        AdditionalPlot, ModVal = self.MakeModFuncs(fluid, Mask)
        v = fluid.u.copy()
        #v[...,0] = ModVal(v[...,0])
        #v[...,1] = ModVal(v[...,1])
        fluid.PlotVectorField(v)
        AdditionalPlot()
    
    def FieldPlot(self, fluid, field, Mask=True, Colorbar=True, colorbar_kwargs={}):
        AdditionalPlot, ModVal = self.MakeModFuncs(fluid, Mask)
        fluid.Plot2D(field=field, additional=AdditionalPlot, mod=ModVal, Colorbar=Colorbar, colorbar_kwargs=colorbar_kwargs)
        
    def StressPlot(self, fluid, Mask = True):			
        AdditionalPlot, ModVal = self.MakeModFuncs(fluid, Mask)
        fluid.StressPlot(f = AdditionalPlot, mod = ModVal)

    def StressDebug(self, fluid, Mask = True):
        AdditionalPlot, ModVal = self.MakeModFuncs(fluid, Mask)
        fluid.StressDebug(f = AdditionalPlot, mod = ModVal)