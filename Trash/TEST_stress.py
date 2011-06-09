


# todo
# finish checking symmetry
# input fake velocity into stress evolution
# setup implicit



#run 'c:\diss her taste, son\peristaltic.py'

import sys, os
root = sys.path[0]
sys.path.append(os.path.join(root, "IB_c\Release"))
import IB_c

import numpy as numpy

from Lookup import *
#from FluidSolve import *
from OB_Fluid import *
from Delta import *

##from Plot import *
import pylab as m
##m.show()
from MatrixUtility import *
from copy import deepcopy
from time import *

from PeristalticPumpGeometry import *

def SecondsToTime(seconds):
    days = int(seconds / (60*60*24))
    seconds -= days * 60*60*24
    hours = int(seconds / (60*60))
    seconds -= hours * 60*60
    minutes = int(seconds / (60))
    seconds -= minutes * 60

    return days, hours, minutes, seconds

# Simulation parameters
DeltaType = 0

N, s, T = 32, 10e5, .25
h = 1. / N
dt =  1000*   10 * h / s**.5
print "dt =", dt

# Setup for saving the simulation output
SaveName = 'Peristaltic_Explicit'
SaveDataFlag, SaveDir, SavePeriod = True, SaveName, 10 #int((T / dt) / 120)
if SaveDataFlag:
    import SaveData as SaveData
    Save = SaveData.SaveDataControl(root, SaveDir)

    
# Peristaltic pump
X = PeristalticPump(Nb_x = N, s = s, Chi = .1)

# Initialize the fluid solver
u = OB_Fluid(dt = dt, N = [N,N,N], dims = [1., 1., 1.], W = 1, B = 1., DeltaType = DeltaType)
#u = Fluid(dt = dt, N = [N,N,N], dims = [1., 1., 1.], DeltaType = DeltaType)




u.u = zeros((N,N,N,3),float64)
#u.u[5,10,3] = [0,0,1]
u.u[:,0:2,0]=[1,0,0]
#u.u[0:2,:,0]=[0,1,0]
#u.u[0,0:2,:]=[0,0,1]
for i in range(N):
    #u.u[:,i,:]=[i*h,0,0]
    u.u[:,:,i]=[0,i*h,0]
for i in range(100):
    Derivatives.TotalDifference_CentralPeriodic (u.N, u.h, u.S, u.dS)
    Derivatives.TotalDifference_CentralPeriodic (u.N, u.h, u.u, u.du)
    u.UpdateS()
    print u.S[5,N/2+5,2]

########
########
########Total_Time = 0
########Times = []
########
########CurStep = 0
########t = 0
########while t < T:
########    print "########################################################################"    
########    print "CurStep =", CurStep, ", t =", t
########    print "u[N/2,N/2,N/4] =", u.u[N/2,N/2,N/4]
########    print "u[N/2,N/2,3N/4] =", u.u[N/2,N/2,3*N/4]    
########    print "u.z.max() =", u.u[:,:,:,2].max(), ", u.z.min() =", u.u[:,:,:,2].min()
########    print
########    print
########
########
########    if abs(X.X).max() > 1:
########        print "Simulation has exploded..."
########        raw_input("")
########        raise
########
########    # Update peristaltic pump tether points
########    X.PeristalticSet(t)
########
########    clock_t = time()
########    X.CalcForceDistribution()
########    FiberToGrid (u.N, u.h, X.Nb, X.hb2, X.X, X.F, u.f, DeltaType)
########    #u.f += Current * flow_source
########    u.FluidSolve(u.f)
########    GridToFiber (u.N, u.h, X.Nb, X.hb2, X.X, u.u, X.U, DeltaType)
########
########    X.X += dt * X.U;
########
########    CurStep += 1
########    t += dt
########
########    Total_Time += time() - clock_t
########    Times.append(Total_Time)
########    
########    seconds = ((T - t) / dt) * (Total_Time / CurStep)
########    days, hours, minutes, seconds = SecondsToTime(seconds)
########    total_seconds = (T / dt) * (Total_Time / CurStep)
########    total_days, total_hours, total_minutes, total_seconds = SecondsToTime(total_seconds)
########
########    print "Average cost per timestep:", Total_Time / CurStep
########    print "Estimated time left:", days, "days,", hours, "hours,", minutes, "minutes,", seconds, "seconds"
########    print "Estimated total time:", total_days, "days,", total_hours, "hours,", total_minutes, "minutes,", total_seconds, "seconds"
########    print "                      (", (Total_Time / CurStep) * (T / dt) / 3600., ") hours, dt =", dt
########    print
########    print
########    print
########
########    
########
########    #if True:
########    #if CurStep % 10 == 0:
########    if SaveDataFlag and CurStep % SavePeriod == 0:
########        m.clf()
########        #z = u.u[:,:,N/2,0]
########        z = u.S[:,:,N/2,0]
########        z = u.S[:,N/2,:,4] - u.S[:,:,N/2,4]
########        print "!!!!!!!!!!!!!!!! ------> ", abs(z).max()
########        m.imshow(z.T, interpolation='bilinear', extent=[0,1,0,1])
########           #vmin = -2, vmax = 10
########        print "%%% ", u.u[:,:,N/2,0].max(), u.u[:,:,N/2,0].min()
########        X.PlotSlice()
########        m.xlim(0, 1)
########        m.ylim(0, 1)
########        #m.show()
########        #raw_input("")
########        #if CurStep == 2: raise
########        Save.SaveStep(u, X, CurStep, SaveFigure = True, Time = Times)
########    
########print
########print "Total simulation time:", Total_Time
