
import numpy as np
import math
 
import matplotlib.pyplot as plt
import sys
import csv

datafile = np.genfromtxt('n_k_JC.OUT')  


theta = (60.0 * np.pi)/180.0
 
n_tilde =  np.zeros(len(datafile),dtype=complex)
for i in range(len(datafile)):
	
	n_tilde[i] = complex(i, i**2)
	print (n_tilde[i])

normal_incidence = ((datafile[:, 1] - 1.0)**2 + datafile[:, 2]**2) /((datafile[:, 1] + 1.0)**2 + datafile[:, 2]**2)
p_polarized = datafile[:, 1]
#for i in range(len(datafile)):   
	#print (normal_incidence[i], '\t', p_polarized[i])         

plt.plot(datafile[:,0], normal_incidence, 'bo') 
plt.show()   	     
 
np.savetxt('expReflectivity.OUT', np.column_stack([datafile[:,0], normal_incidence]), delimiter='\t') 








