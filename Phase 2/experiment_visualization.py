from Utils.define_experiment import experiment1
import subfunctions as sf
from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt


experiment, end_event = experiment1()
alpha_dist = experiment['alpha_dist']
alpha_deg = experiment['alpha_deg']



alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate') #fit the cubic spline


x = np.linspace(alpha_dist[0], alpha_dist[-1], 100)


y= []
for i in range(len(x)):
    y.append(alpha_fun(x[i]))


plt.plot(x, y,label = 'Cubic Spline')
plt.plot(alpha_dist, alpha_deg, '*', label = 'Data Points', )
plt.xlabel('Distance (m)')
plt.ylabel('Alpha (degrees)')
plt.title('Alpha vs Distance')
plt.legend()
plt.show()