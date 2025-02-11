import numpy as np
import math

data = np.genfromtxt("data/Obergfell_transmissivity.csv",delimiter=",")

width = data[0,:].size
height = data[:,0].size
print(width-1)

x = []
y = []
z = []
for i in range(height-1):
    x.extend([data[i+1,0] for j in range(width-1)])
    y.extend(data[0,1:])
    z.extend([data[i+1,j+1] for j in range(width-1)])


np.savetxt('data/Obergfell_transmissivity_gle.OUT', np.column_stack([x,y,z]), delimiter='\t')
