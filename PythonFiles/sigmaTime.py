
import numpy as np
import math
 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import csv

print("\nStart Ndione's plotting script")

exp =  np.genfromtxt("/home/ndione/Desktop/goldflash/expDataReanalyzed/750nmsigma.txt", comments="#")

x_exp = exp[:, 0]
x_exp_errorbar = exp[:, 1]
y_exp = exp[:, 2]
y_exp_errorbar = exp[:, 3]

DFT = np.genfromtxt("sigma_r_1_65.OUT", comments="#")
x_DFT = DFT[:, 0] 
y_DFT = DFT[:, 1]

convert = 1e-15
Ndione_0_8 = np.genfromtxt("/home/ndione/NetBeansProjects/DF_Program/var/NewXFEL/newDOS/0_8MJ_kg/1_65eV/sigma_const_ph.OUT", comments="#")
x_Ndione_0_8 = Ndione_0_8[:, 0] 
y_Ndione_0_8 = Ndione_0_8[:, 1] * convert
Ndione_1_02 = np.genfromtxt("/home/ndione/NetBeansProjects/DF_Program/var/NewXFEL/newDOS/1_02MJ_kg/1_65eV/sigma_const_ph.OUT", comments="#")
x_Ndione1_02  = Ndione_1_02[:, 0]
y_Ndione_1_02 = Ndione_1_02[:, 1] * convert
Ndione_1_2 = np.genfromtxt("/home/ndione/NetBeansProjects/DF_Program/var/NewXFEL/newDOS/1_2MJ_kg/1_65eV/sigma_const_ph.OUT", comments="#")
x_Ndione_1_2 = Ndione_1_2[:, 0]
y_Ndione_1_2 = Ndione_1_2[:, 1] * convert
Ndione_1_6 = np.genfromtxt("/home/ndione/NetBeansProjects/DF_Program/var/NewXFEL/newDOS/1_6MJ_kg/1_65eV/sigma_const_ph.OUT", comments="#")
x_Ndione_1_6 = Ndione_1_6[:, 0]
y_Ndione_1_6 = Ndione_1_6[:, 1] * convert


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}')
fig, axis = plt.subplots(figsize=(4,3))
axis.set_xlabel('time ($\si{\\femto\second}$)')
axis.set_ylabel('$\sigma_1$ ($\SI{1e15}{\second^{-1}}$)')

ndione_plot = axis.plot(x_Ndione_1_6, y_Ndione_1_6, color='red', linestyle='-', label='Present work')
DFT_plot = axis.plot(x_DFT, y_DFT, color='green', linestyle='--', label='DFT')
exp_plot = axis.errorbar(x_exp, y_exp, xerr = x_exp_errorbar,  yerr = y_exp_errorbar,   color='b', fmt='o', label='exp.', capsize = 1.2)


axis.set_xlim(left=-500.0, right=2200.0)
#axis.set_ylim(0.0, 10000.0*scale)
axis.set_ylim(0.0, 6.0)
axis.tick_params(axis='y')

axis.legend()
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

print("\nEnd Ndione's plotting script")
