from typing import NamedTuple

import numpy as np
from pcsaft import *
from param import *

method = PCSAFT()
param = Param
param.m = np.array([1, 2, 3, 4])
param.e = np.array([200, 300, 100, 180])
param.s = np.array([2, 3, 4, 5])
param.k_ij = np.zeros([4, 4])
comp = np.array([0.3, 0.2, 0.1, 0.4])
T = 298.15
p = 101325

#ztest = method.compute_z(0.04, T, comp, param)
rhotest=method.compute_density(T,p,comp,0,param)

print(rhotest)
