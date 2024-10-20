import numpy as np
from matplotlib import pyplot as plt
from define_experiment import experiment1
import subfunctions as sf
import pandas as pd


# Load the rover dictionary from subfunctions
rover = {
  'wheel_assembly': {
      'wheel': {'radius': 0.30, 'mass': 1.0
                },
      
      'motor': {'torque_stall': 170, 'torque_noload': 0, 'speed_noload': 3.80, 'mass': 5.0, 
                'effcy_tau' : np.array([0, 10, 20, 40, 70, 165]) ,
     'effcy' : np.array([0, .55, .75, .71, .50, .05])
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

# Load experiment and end_event data
experiment, end_event = experiment1()

# Update end_event parameters as specified
end_event['max_distance'] = 1000
end_event['max_time'] = 10000
end_event['min_velocity'] = 0.01

# Run simulation
rover = sf.simulate_rover(rover, planet, experiment, end_event)

# Create figure with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

# Position vs time plot
ax1.plot(rover['telemetry']['Time'], rover['telemetry']['position'])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (m)')
ax1.set_title('Position vs Time')
ax1.grid(True)

# Velocity vs time plot
ax2.plot(rover['telemetry']['Time'], rover['telemetry']['velocity'])
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Velocity (m/s)')
ax2.set_title('Velocity vs Time')
ax2.grid(True)

# Power vs time plot
ax3.plot(rover['telemetry']['Time'], rover['telemetry']['power'])
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Power (W)')
ax3.set_title('Power vs Time')
ax3.grid(True)

plt.tight_layout()
plt.show()

# Print telemetry data
#print("\nRover Telemetry Data:")
#print(f"Completion Time: {rover['telemetry']['completion_time']:.2f} s")
#print(f"Distance Traveled: {rover['telemetry']['distance_traveled']:.2f} m")
#print(f"Maximum Velocity: {rover['telemetry']['max_velocity']:.2f} m/s")
#print(f"Average Velocity: {rover['telemetry']['average_velocity']:.2f} m/s")
#print(f"Battery Energy Used: {rover['telemetry']['battery_energy']:.2f} J")
#print(f"Energy per Distance: {rover['telemetry']['energy_per_distance']:.2f} J/m")

##print(rover['telemetry'])
