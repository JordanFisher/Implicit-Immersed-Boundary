from numpy import array

HoldPythonArrayType = type([])
def EnsureNumpyArray(l):
    """Convert l to a Numpy array if necessary, then return."""

    if type(l) == HoldPythonArrayType:
        return array(l)
    else:
        return l

def cross(v, w):
    """Return the cross product v X w."""

    z = v.copy()
    z[0] = v[1] * w[2] - v[2] * w[1]
    z[1] = -v[0] * w[2] + v[2] * w[0]
    z[2] = v[0] * w[1] - v[1] * w[0]

    return z
    
