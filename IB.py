import os, sys

from os.path import join, dirname

path = dirname(__file__)

sys.path.append(join(path, 'Fiber'))
sys.path.append(join(path, 'Fluid'))
sys.path.append(join(path, 'IB_Utility'))
sys.path.append(join(path, 'Implicit_Utility'))
sys.path.append(join(path, 'Lookups'))
sys.path.append(join(path, 'Sim_Utility'))
sys.path.append(join(path, 'Plot_Utility'))
sys.path.append(join(path, 'FastEval'))
sys.path.append(join(join(join(path, 'FastEval'), 'IB_c'), 'Release'))
sys.path.append(join(join(join(path, 'FastEval'), 'Panel_c'), 'Release'))