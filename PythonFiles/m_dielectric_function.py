import numpy as np
import math

data = np.genfromtxt("data/JC_300K_CU_n,k.csv", delimiter=',')

x = data[:,0]
epsilon = [complex(data[i,1],data[i,2])**2 for i in range(x.size)]
y1 = [epsilon[i].real for i in range(x.size)]
y2 = [epsilon[i].imag for i in range(x.size)]

np.savetxt('data/JC_300K_CU_epsilon.OUT', np.column_stack([ x, y1, y2]), delimiter='\t')
