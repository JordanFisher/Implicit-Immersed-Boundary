import numpy as numpy

import Panel
import Decomposition as D
import DecompGroup as G

from numpy import *

from os import path

def DecompFolderRoot(root, fluid):
    return path.join(root, fluid.Params_ToString())

def InitDecomposition(Decomp, file_root, fluid, PanelWidth, Terms):
    SomeFilesExist = False
    
    Decomp.Initialize(Terms, PanelWidth)    
    file_root = DecompFolderRoot(file_root, fluid)
    file_root = path.join(file_root, 'Panel_Width_' + str(PanelWidth))

    for l in range(Terms):
        _N = fluid.N[0] / 2 + PanelWidth / 2 + 2
        u = zeros((3,3,_N,_N,_N))

        name = path.join(file_root, str(l) + '_00.npy')
        try:
            loaded = numpy.load(name)
        except:
            if SomeFilesExist:
                raise Exception("Not enough terms exist.")
            else:
                raise Exception("Decomposition doesn't exist.")
        
        Decomp.Set_AB_xx(l, loaded)

        name = path.join(file_root, str(l) + '_01.npy')
        Decomp.Set_AB_xy(l, numpy.load(name))

        SomeFilesExist = True

def DecompositionGroup(file_root, fluid, Terms, MinPanelWidth):
    """Return a Decomposition Group loaded from the folder located at file_root."""
    
    SomeLevelsExist = False # A flag set to True once we know at least one level of the decomposition group exists
    
    N = fluid.N[0]

    # Create a list to hold the different levels of the decomposition group. The first two levels are empty.
    UField = [None, None] 

    PanelWidth = N / 4
    while PanelWidth > MinPanelWidth:
        Decomp = D.Decomposition()

        try:
            InitDecomposition(Decomp, file_root, fluid, PanelWidth, Terms)
        except Exception, (msg):
            if SomeLevelsExist:
                raise Exception("Not enough levels exist.")
            else:
                raise Exception(msg)

        SomeLevelsExist = True
        
        UField.append(Decomp)
        PanelWidth /= 2

    Group = G.DecompGroup()

    Group.Initialize(len(UField), N)
    for i in range(len(UField)):
        Group.Set_UField(i, UField[i])

    return Group

