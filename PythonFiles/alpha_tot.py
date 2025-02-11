
import numpy as np
import math
 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import csv

datei = np.genfromtxt('n_k_JC.OUT')

data = np.genfromtxt('n_k_Werner.OUT')  

light_speed = 299792458

eV_to_Hz =  2.417989262e14

hbar_eV = 6.582119569e-16

omega_JC = datei[:,0]  / hbar_eV

omega_Werner = data[:, 0]  / hbar_eV

k_JC = datei[:, 2]

k_Werner = data[:, 2]

alpha_JC = (2.0 * omega_JC * k_JC)/light_speed

alpha_Werner = (2.0 * omega_Werner * k_Werner) / light_speed
 

for i in range(len(datei)):

		print ("JC = ", alpha_JC[i], '\t', "Werner = ", alpha_Werner[i])         

np.savetxt('alpha_tot_JC.OUT', np.column_stack([datei[:,0], alpha_JC]), delimiter='\t')   

np.savetxt('alpha_tot_Werner.OUT', np.column_stack([data[:,0], alpha_Werner]), delimiter='\t')  

figure1 = plt.figure()
plt.xlabel('Energy [eV]')
plt.ylabel(r'$\alpha_{tot}$ [$\times 10^{7}\,m^{-1}$]')
plt.plot(datei[:,0], alpha_JC*1e-7, 'go',  markersize=3)
figure1.savefig('JC.pdf')


figure2 = plt.figure()
plt.xlabel('Energy [eV]')
plt.ylabel(r'$\alpha_{tot}$ [$\times 10^{7}\,m^{-1}$]')
plt.plot(data[:,0], alpha_Werner*1e-7, 'ro',  markersize=3)
figure2.savefig('Werner.pdf')
plt.show()  

print ("\nScript run correctly!")







