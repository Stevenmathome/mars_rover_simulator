import numpy as np
import matplotlib.pyplot as plt
from subfunctions import F_net, get_mass, get_gear_ratio

# Dictionaries
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

# Generate rolling resistance coefficients
Crr_array = np.linspace(0.01, 0.4, 25)

def find_terminal_velocity(Crr):
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    wheel_radius = rover['wheel_assembly']['wheel']['radius']
    speed_noload = rover['wheel_assembly']['motor']['speed_noload']
    torque_stall = rover['wheel_assembly']['motor']['torque_stall']
    
    max_velocity = speed_noload * wheel_radius / gear_ratio
    
    # Create an array of velocities
    velocities = np.linspace(0, max_velocity, 1000)
    
    # Calculate motor angular velocity
    motor_omega = velocities * gear_ratio / wheel_radius
    
    # Calculate motor torque (assuming linear torque-speed curve)
    motor_torque = torque_stall * (1 - motor_omega / speed_noload)
    
    # Calculate driving force
    driving_force = motor_torque * gear_ratio / wheel_radius
    
    # Calculate resistive force (rolling resistance)
    mass = get_mass(rover)
    resistive_force = Crr * mass * planet['g']
    
    # Find where driving force equals resistive force
    force_diff = driving_force - resistive_force
    terminal_velocity_index = np.argmin(np.abs(force_diff))
    terminal_velocity = velocities[terminal_velocity_index]
    
    print(f"Crr: {Crr}, Terminal velocity: {terminal_velocity}")
    
    return terminal_velocity

# Calculate maximum velocity for each rolling resistance coefficient
v_max = np.array([find_terminal_velocity(Crr) for Crr in Crr_array])

# Print v_max to check its contents
print("v_max array:", v_max)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(Crr_array, v_max)
plt.xlabel('Coefficient of Rolling Resistance')
plt.ylabel('Maximum Velocity (m/s)')
plt.title('Maximum Rover Velocity vs. Rolling Resistance Coefficient')
plt.grid(True)
plt.show()