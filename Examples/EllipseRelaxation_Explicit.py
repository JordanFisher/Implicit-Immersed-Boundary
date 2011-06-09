"""A simple example of the explicit Immersed Boundary Method.
An ellipse is constructed and allowed to relax to a circular configuration.
The explicit method used is FE/BE (Forward Euler/Backward Euler)
"""

import IB
import pylab as m
from numpy import array, zeros, ones, float64

from Fluid import Fluid, PlotStyle
from Ellipse import Ellipse
from SimTime import SimTime
import SaveData as SaveData

# Simulation parameters
N, s, T = 64, 10e6, .102
h = 1. / N

dt = .00000003

Time = SimTime(dt = dt, T = T)

# Setup for saving the simulation output
SavePeriod = 100 #int((T / dt) / 120)
Save = SaveData.SaveDataControl(('EllipseRelaxation', 'Explicit'))

    
# Ellipse
X = Ellipse(2 * N, s, [.5,.5], [.8,.5], [.5,.6])
X.UseNativePython = False # If true native python code is used instead of C code

# Initialize the fluid solver
u = Fluid(dt = dt, N = [N,N], dims = [1., 1.])

# The simulation loop        
Time.StartTimer()
while Time.t < Time.T:
    Time.PrintStepInfo()

    X.CheckForExplosion()

    # The FE/BE method
    # First calculate the fiber force F (stored in X.F)
    X.CalcForceDistribution()

    # Spread F to the Eulerian grid, storing the result in the force field u.f
    X.ToGrid(u)

    # Evolve the fluid forward one time step
    u.UpdateFluid(u.f, CalcPressure = True)

    # Interpolate the fluid velocity to the fiber (stored in X.U)
    X.FromGrid(u)

    # Update the fiber position
    X.X += dt * X.U

    # Update the time
    Time.Incr()
    Time.PrintTimeInfo()

    # Plot the fluid/fiber and save to file
    if Time.CurStep % SavePeriod == 0:
        m.clf() # Clear the plotter
        u.Plot2D(PlotStyle.pressure) # Plot the fluid pressure
        X.Plot2D() # Plot the ellipse
        
        Save.SaveStep(u, X, Time.CurStep, SaveFigure = True)
