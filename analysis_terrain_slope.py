import numpy as np
from scipy.optimize import brentq, minimize_scalar
from subfunctions import F_net
import matplotlib.pyplot as plt

def v_max_at_slope(slope, rover, planet, Crr):
    def f(omega):
        return F_net(omega, slope, rover, planet, Crr)
    
    omega_nl = rover['wheel_assembly']['motor']['speed_noload']
    
    try:
        # First, try brentq
        omega_max = brentq(f, 0, omega_nl)
    except ValueError:
        # If brentq fails, try to find the maximum velocity
        def neg_v(omega):
            return -omega * rover['wheel_assembly']['wheel']['radius'] / rover['wheel_assembly']['speed_reducer']['diam_gear'] * rover['wheel_assembly']['speed_reducer']['diam_pinion']
        
        result = minimize_scalar(neg_v, bounds=(0, omega_nl), method='bounded')
        omega_max = result.x
    
    return omega_max * rover['wheel_assembly']['wheel']['radius'] / rover['wheel_assembly']['speed_reducer']['diam_gear'] * rover['wheel_assembly']['speed_reducer']['diam_pinion']

rover = {
    'wheel_assembly': {
        'wheel': {'radius': 0.30, 'mass': 1.0},
        'speed_reducer': {'type': 'reverted', 'diam_pinion': 0.04, 'diam_gear': 0.07, 'mass': 1.5},
        'motor': {'torque_stall': 170, 'torque_noload': 0, 'speed_noload': 3.80, 'mass': 5.0}
    },
    'chassis': {'mass': 659},
    'science_payload': {'mass': 75},
    'power_subsys': {'mass': 90}
}

planet = {'g': 3.72}

slope_array_deg = np.linspace(-10, 35, 25)
Crr = 0.2

v_max = np.array([v_max_at_slope(slope, rover, planet, Crr) for slope in slope_array_deg])

plt.figure(figsize=(10, 6))
plt.plot(slope_array_deg, v_max)
plt.xlabel('Terrain Slope (degrees)')
plt.ylabel('Maximum Velocity (m/s)')
plt.title('Maximum Rover Velocity vs. Terrain Slope')
plt.grid(True)
plt.show()