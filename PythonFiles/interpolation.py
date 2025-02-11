import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import interpolate
import sys
import io



data =  np.genfromtxt('alpha.OUT')


func1 = interp1d(data[:, 0], data[:, 1], kind='linear')

func2 = interp1d(data[:, 0], data[:, 1], kind='cubic')

x_interp = np.linspace(0.5, 6, num=100, endpoint=True) 
 

#interpolate.spl(x, y, s): perform a spline interpolation
#Argument s is used to specify the amount of smoothing to perform during the spline fit. 
#The default value of s if not provide is s=n-sqrt(n), n being the number of data points to fit

##ynew = interpolate.splrep(data[:, 0], data[:, 1]) 
##func3 =  interpolate.splev(xnew, ynew, der=0)


plt.plot(data[:, 0], data[:, 1], 'o', x_interp, func, '<', xnew, func2(xnew))
plt.legend(['data', 'linear', 'cubic'], loc='best')
plt.show()
