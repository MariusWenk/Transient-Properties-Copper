
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from unit import *

data = np.genfromtxt('dynamics.OUT')

time = data[:, 0]

n_sp = data[:, 1]

n_d = data[:, 2]

Tph = data[:, 3]

wp =  (np.sqrt((n_sp*e_charge*e_charge)/(m_e*eps0))/(2.0*np.pi)) *Hz_to_eV

wp_square = wp*wp

nu_eph = (0.129e15*(Tph/300.))*Hz_to_eV

nu_ee = 0.36*n_d*unit_volume*((n_d[0]-n_d)*unit_volume)*1e15*Hz_to_eV*0.4

nu_tot = nu_eph + nu_ee

eps_r = 3.104-(wp_square/(omega_probe**2+nu_tot**2))

eps_i = (nu_tot*wp_square)/(omega_probe*(omega_probe**2+nu_tot**2))

eps = np.zeros(len(data), dtype=np.complex)

refrac_index  = np.zeros(len(data), dtype=np.complex)

n_reflectivity  = np.zeros(len(data), dtype=float)

s_reflectivity  = np.zeros(len(data), dtype=float)

p_reflectivity  = np.zeros(len(data), dtype=float)

for idx in range(len(data)):

	eps[idx] = np.complex(eps_r[idx], eps_i[idx])

	refrac_index[idx] = np.sqrt(eps[idx])

	n_reflectivity[idx] = np.abs((1.0-refrac_index[idx])/(1.0+refrac_index[idx]))**2
 
for idx in range(len(data)):

	p_reflectivity[idx] = abs((np.cos(theta)*refrac_index[idx]**2 - np.sqrt(refrac_index[idx]**2-np.sin(theta)**2)) / (np.cos(theta)*refrac_index[idx]**2 + np.sqrt(refrac_index[idx]**2-np.sin(theta)**2)))**2
	#print (p_reflectivity[idx])

for idx in range(len(data)):

	s_reflectivity[idx] = abs(0)**2


np.savetxt("radfile.OUT", np.column_stack([time, n_reflectivity]), delimiter = "\t")
plt.xlim(-500, 6000)
plt.plot(time, p_reflectivity)
#plt.show()
 







