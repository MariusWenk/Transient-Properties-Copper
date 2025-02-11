import numpy as np
import math

Hartree_to_eV = 27.2114

Emax = 4.3114e+01
#Emax = 10

numbers = np.linspace(-10,53,10000)
#bsp = 7*(3/2)/((Emax+9.361)**(3/2))
#bd = 10*(3/2)/((-2.1+9.361)**(3/2))
bsp = 1.4252/Hartree_to_eV
bd = 20.8712/Hartree_to_eV
print(bsp)
print(bd)

sp_DOS_stretching = 0.7573
d_DOS_stretching = 1.0372

sp_data = []
d_data = []
for i in numbers:
    if i <= Emax and i>= -9.361:
        sp_data.append(bsp*np.sqrt(i+9.361))
    else:
        sp_data.append(0)
        
    if i <= -2.1 and i>= -9.361:
        d_data.append(bd*np.sqrt(i+9.361))
    else:
        d_data.append(0)
        
sp_data = np.array(sp_data)
d_data = np.array(d_data)
total_data = np.array(sp_data + d_data)
 

np.savetxt('data/DOS_easy_Cu.in', np.column_stack([ numbers, sp_data, d_data, total_data]), delimiter='\t')
np.savetxt('../input/DOS_Cu/DOS_easy_Cu.in', np.column_stack([ numbers/Hartree_to_eV, sp_data*Hartree_to_eV/sp_DOS_stretching, d_data*Hartree_to_eV/d_DOS_stretching, total_data]), delimiter='\t')
