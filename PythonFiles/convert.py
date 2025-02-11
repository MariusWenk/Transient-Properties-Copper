import numpy as np
import math

data = np.genfromtxt("../input/coupling_SI/Medvedev_coupling_SI.in")

K_in_au = 3.16681542254311e-6
Joule_to_Hartree = 1.0/4.35974417e-18
s_to_hbar_E_H = 4.1341e16
m_to_a0 = 1.8898e10
eph_in_au = (Joule_to_Hartree)/(s_to_hbar_E_H*K_in_au*m_to_a0** 3.0)
x = data[:,0]*K_in_au 
y = data[:, 1]*eph_in_au
 

np.savetxt('../input/coupling_au/Medvedev_coupling_au.in', np.column_stack([ x, y]), delimiter='\t')
