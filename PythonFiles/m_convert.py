import numpy as np
import math

data = np.genfromtxt("../input/DOS_Cu/Cu_l_proj_dos.dat")

Hartree_to_eV = 27.2114

x = data[:,0] / Hartree_to_eV
y1 = (data[:, 2]+data[:, 3]) * Hartree_to_eV
y2 = data[:, 4] * Hartree_to_eV
#y3 = data[:, 1]
y3 = y1+y2
 

np.savetxt('../input/DOS_Cu/Cu_l_proj_dos.in', np.column_stack([ x, y1, y2, y3]), delimiter='\t')
