import Xenon as Xe # Loads the polynomials in a python script. 
import numpy as np
import matplotlib.pyplot as plt

def total_response(mass, c1):
    Xe.PP['mdm']=mass
    Xe.PP['c1']=c1
    return sum(Xe.ret_ipol(Xe.PP)['/DM/XENON/counts'])

def find_coeff(mass):
    logmass = np.log10(mass)
    #print(logmass)
    c1 = -5.0
    counts=total_response(logmass, c1)
    return 1e-5*np.sqrt(2.3/counts)

mass_vals = np.logspace(1.04,4.0,1000)
coeffs = []
for i in range(len(mass_vals)):
    coeffs.append(find_coeff(mass_vals[i]))
plt.loglog(mass_vals, coeffs)
data_30=np.genfromtxt('exclusion_Xe_30keV.dat', delimiter=' ')
plt.loglog(data_30[:,0], data_30[:,1])

plt.show()
