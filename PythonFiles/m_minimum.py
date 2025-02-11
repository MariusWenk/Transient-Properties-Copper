import numpy as np
import math

data = np.genfromtxt("data/FullDynamics2_filtered.OUT")

Hartree_to_eV = 27.2114

x = data[:,0]
y = data[:,1]

ymin = np.amin(y)
index = np.where(y == ymin)
print (f"Minimum at {x[index][0]} fs, {(y[0]-ymin)*1e-6} 1/cm^3 below initial value")
#print (np.where(y == y.min()))

