
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import csv

argon = np.genfromtxt('JC_Au_refractive_index_n_k.dat')  


factor1 = (1.0 - argon[:,1])**2 + (argon[:, 2])**2
factor2 = (1.0 + argon[:, 1])**2 + (argon[:, 2])**2
reflectivity = factor1/factor2

transmission  = 1.0 - reflectivity

real_diel = (argon[:, 1])**2 - (argon[:, 2])**2
imag_diel = 2.0 * argon[:, 1] * argon[:, 2]
 

for i in range(len(argon)):   
	print transmission[i]         

plt.plot(argon[:,0], transmission, 'bo') 
plt.show()   	     
#plt.xlim(1.0, 6.0)
#plt.ylim(0.0, 7.0)
np.savetxt('JC_Au_transmission.dat', np.column_stack([argon[:,0], transmission]), delimiter='\t') # save your results! :) 








