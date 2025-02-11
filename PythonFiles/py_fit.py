import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

def objectiveFunction(x, a, b, c, d):
	return a*x**3+b*x**2+c*x+d

data = np.genfromtxt("f_below_0_8")
xdata =data[:, 0]
ydata=data[:, 1]

popt, pcov = curve_fit(objectiveFunction, xdata, ydata)

params=['a', 'b', 'c', 'd']

for i in range(len(popt)):
	print(params[i], "=", popt[i])

print ("\npcov ", pcov)

plt.plot(xdata, ydata, 'r-', label='orignial data',)
plt.plot(xdata, objectiveFunction(xdata, *popt), 'b--', label='fit')
plt.legend()
plt.xlabel("temperature, K")
plt.ylabel("energy, eV")
plt.show()
