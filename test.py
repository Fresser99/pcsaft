import matplotlib.pyplot as plt
import numpy as np
from pcsaft import *

from param import *

param = Param
param.CAS = np.array(['seg-IB-R', 'seg-IP-R'])
param.type = np.array([1, 1])
T = np.linspace(179, 373, 100)
copolyfrac = np.array([0.98, 0.02])
method = PCSAFT()
den = np.array([])
cp=np.array([])
lbd=np.array([])
for t in T:
    den = np.append(den, method.compute_solid_density_Askadskii_Matveev(param, t, copolyfrac))
    cp=np.append(cp,method.compute_solid_heat_capacity_Bicerano(param,t,copolyfrac))
    lbd=np.append(lbd,method.compute_solid_thermal_conductivity_Askadskii_Matveev(param,copolyfrac,t))
fig = plt.figure()
#plt.plot(T, den)
plt.plot(T,lbd)
plt.show()
