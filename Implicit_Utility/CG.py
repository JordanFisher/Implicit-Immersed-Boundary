import numpy as numpy

def Dot(a, b):
    return (a*b).sum()

def CG (x, b, f, Tolerance = 10e-8, MaxReps = 1000, ShowProgress = False, Callback = None):
    """ Conjugate Gradient:
        Attempts to solve the system f(y) = b for y
        x is the initial guess"""

    r = b - f(x)

    if abs(r).max() == 0:
        return x, 0
    
    w = -r;

    z = f(w)

    a = Dot(r, w) / Dot(w, z)
    
    print "!!!", abs(x).max(), a, abs(w).max()
    x = x + a * w
    print "!!!", abs(x).max()

    Reps = 1

    for q in range(MaxReps):
        r = r - a * z

        if ShowProgress:
            print "Step", q, ", ||r|| =", Dot(r,r)**.5, ", max(|r|) =", abs(r).max()
            if Callback != None:
                Callback(x, r)

        if abs(r).max() < Tolerance:
            break
##        if Dot(r,r)**.5 < Tolerance:
##            return x

        B = Dot(r, z) / Dot(w, z)

        w = -r + B * w

        z = f(w)

        a = Dot(r, w) / Dot(w, z)

        x = x + a * w

        Reps += 1

    return x, Reps
