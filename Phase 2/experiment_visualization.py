from define_experiment import experiment1
import subfunctions as sf
from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt

# define experiment and end_event along with the angle at different distances
experiment, end_event = experiment1()
alpha_dist = experiment['alpha_dist']
alpha_deg = experiment['alpha_deg']


#create a interpolation function to interpolate the data
alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate') #fit the cubic spline

# create a continuous range of distance values for plotting
x = np.linspace(alpha_dist[0], alpha_dist[-1], 100)

#get the corresponding alpha values
y= []
for i in range(len(x)):
    y.append(alpha_fun(x[i]))

#plot alpha vs distance as a continous function and data points
plt.plot(x, y,label = 'Cubic Spline')
plt.plot(alpha_dist, alpha_deg, '*', label = 'Data Points', )
plt.xlabel('Distance (m)')
plt.ylabel('Alpha (degrees)')
plt.title('Alpha vs Distance')
plt.legend()
plt.show()