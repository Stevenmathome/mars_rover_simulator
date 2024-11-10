from define_experiment import experiment1
import subfunctions as sf
from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt


rover = {
  'wheel_assembly': {
      'wheel': {'radius': 0.30, 'mass': 1.0
                },
      
      'motor': {'torque_stall': 170, 'torque_noload': 0, 'speed_noload': 3.80, 'mass': 5.0, 
                'effcy_tau' : np.array([0, 10, 20, 40, 70, 165]) ,
     'effcy' : np.array([0, .55,.75,.71,.50,.05])
                },
      
      'speed_reducer': {'type': 'reverted', 'diam_pinion': 0.04, 'diam_gear': 0.07, 'mass': 1.5}
  },
  'chassis': {'mass': 659
              },
  
  'science_payload': {'mass': 75
                      },
  
  'power_subsys': {'mass': 90
                   }
}

planet = {'g': 3.72}

# define experiment and end_event
experiment , end_event = experiment1()

# get effcy tay and effcy from rover
effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']   
effcy = rover['wheel_assembly']['motor']['effcy']

#create effcy function to interpolate the data
effcy_fun = interp1d(effcy_tau, effcy, kind = 'cubic', fill_value='extrapolate') #fit the cubic spline
# create a continuous range of torque values for plotting
x = np.linspace(effcy_tau[0], effcy_tau[-1], 100)

#plot effcy vs torque as a continous function and data points
plt.plot(x, effcy_fun(x)*100,label = 'Cubic Spline')
plt.plot(effcy_tau, effcy*100, '*', label = 'Data Points', )
plt.xlabel('Torque (N-m)')
plt.ylabel('Efficiency (%)')
plt.title('Efficiency (%) vs Torque')
plt.legend()
plt.show()   