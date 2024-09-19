from typing import NamedTuple

import numpy as np
from pcsaft import *


class Param(NamedTuple):
    s: np.ndarray
    e: np.ndarray
    m: np.ndarray
    k_ij: np.ndarray


method = PCSAFT()
param=Param
param.m = np.array([2, 3, 1.5, 3])
param.e = np.array([10, 11, 21, 31])
param.s = np.array([200, 171, 300, 111])
param.k_ij = np.zeros([4,4])
comp = np.array([0.8, 0.1, 0.05, 0.05])
T = 298.15
p = 101325
ztest=method.compute_z(0.03,T,comp,param)
print(ztest)
