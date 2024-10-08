# Pat
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import bisect
from subfunctions import F_net, tau_dcmotor, get_mass, get_gear_ratio, F_gravity, F_rolling

# Constants
Crr_array = np.linspace(0.01, 0.4, 25)
slope_array_deg = np.linspace(-10, 35, 25)

# Create meshgrid for the 3D plot
CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)

# Create a matrix to store maximum speeds
VMAX = np.zeros(np.shape(CRR))

# Parameters for the rover and planet
rover = {
    'wheel_assembly': {
        'wheel': {'radius': 0.30, 'mass': 1.0},
        'motor': {'torque_stall': 170, 'torque_noload': 0, 'speed_noload': 3.80, 'mass': 5.0},
        'speed_reducer': {'type': 'reverted', 'diam_pinion': 0.04, 'diam_gear': 0.07, 'mass': 1.5}
    },
    'chassis': {'mass': 659},
    'science_payload': {'mass': 75},
    'power_subsys': {'mass': 90}
}

planet = {'g': 3.72}  

# Loop through each combination of rolling resistance and slope
N = np.shape(CRR)[0]
for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i, j])
        slope_sample = float(SLOPE[i, j])

        # Define a function for the root-finding method (net force equals zero)
        def net_force(omega):
            return F_net(omega, slope_sample, rover, planet, Crr_sample)

        # Use bisection or another method to find the root
        try:
            # Initial guesses based on motor speed range (adjust if needed)
            omega_min = 0  # No movement
            omega_max = rover['wheel_assembly']['motor']['speed_noload']
            
            # Find the speed where the net force is zero (steady-state speed)
            omega_solution = bisect(net_force, omega_min, omega_max)
            VMAX[i, j] = omega_solution * rover['wheel_assembly']['wheel']['radius']  # Convert to translational speed
        except ValueError:
            # In cases where no solution is found, set the value to NaN
            VMAX[i, j] = np.nan

# Plot the surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(CRR, SLOPE, VMAX, cmap='viridis')

# Label axes
ax.set_xlabel('Rolling Resistance Coefficient')
ax.set_ylabel('Terrain Slope (degrees)')
ax.set_zlabel('Max Rover Speed (m/s)')
ax.set_title('Max Speed of Rover Over Various Terrain Conditions')

plt.show()
