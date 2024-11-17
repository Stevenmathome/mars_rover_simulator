"""###########################################################################
# This file is the top-level script for running a simulation of the EDL system.
# It loads the appropriate info into the workspace and calls simulate_edl.
# After returning from simulate_edl, it graphs some of the simulation
# trajectory data.
#
#   Created by: Former Marvin Numerical Methods Engineering Team
#   Last Modified: 28 October 2023
###########################################################################"""

import numpy as np
import matplotlib.pyplot as plt
from define_edl_system import *
from subfunctions_EDL import *
from define_planet import *
from define_mission_events import *
from redefine_edl_system import *

# *************************************
# load dicts that define the EDL system (includes rover), planet,
# mission events and so forth.
edl_system = define_edl_system_1()
mars = define_planet()
mission_events = define_mission_events()


# Overrides what might be in the loaded data to establish our desired
# initial conditions
edl_system['altitude'] = 11000    # [m] initial altitude
edl_system['velocity'] = -590     # [m/s] initial velocity

edl_system['heat_shield']['ejected'] = True

edl_system['sky_crane']['on'] = False

edl_system['parachute']['deployed'] = True   # our parachute is open
edl_system['parachute']['ejected'] = False   # and still attached

edl_system['speed_control']['on'] = False
edl_system['position_control']['on'] = False

edl_system['rover']['on_ground'] = False # the rover has not yet landed
edl_system['parachute']['diameter'] = 18.0

tmax = 2000   # [s] maximum simulated time


points = (19-14)/.5
x = np.linspace(14,19, int(points)+1)

time = []
rover_speed = []
rover_sucess = []
altitude=[]
safe_landing_vecolcity = 1
for i in range(len(x)):
    
    edl_system = redefine_edl_system(edl_system)
    edl_system['parachute']['diameter'] = x[i]
    # the simulation. changing last argument to false turns off message echo
    
    [t, Y, edl_system] = simulate_edl(edl_system, mars, mission_events, tmax, True,Modified=True)
    time.append(t[-1])
    rover_speed.append(Y[0, -1])
    altitude.append(Y[1, -1])
    if abs(rover_speed[-1]) <=abs(safe_landing_vecolcity) and altitude[-1] >= 4.5:
        rover_sucess.append(1)
    else:
        rover_sucess.append(0) 
    

# Create the 3x1 array of plots
fig, axs = plt.subplots(3, 1, figsize=(8, 12))

# Plot 1: Termination time vs. parachute diameter
axs[0].plot(x, time, marker='o')
axs[0].set_title("Simulated Time vs. Parachute Diameter")
axs[0].set_xlabel("Parachute Diameter (m)")
axs[0].set_ylabel("Termination Time (s)")
axs[0].set_xticks(x)
axs[0].grid()

# Plot 2: Rover speed at termination vs. parachute diameter
axs[1].plot(x, rover_speed, marker='o')
axs[1].set_title("Rover Speed at Termination vs. Parachute Diameter")
axs[1].set_xlabel("Parachute Diameter (m)")
axs[1].set_ylabel("Rover Speed (m/s)")
axs[1].set_xticks(x)
axs[1].grid()

# Plot 3: Rover landing success vs. parachute diameter
axs[2].plot(x, rover_sucess, marker='o')
axs[2].set_title("Rover Landing Success vs. Parachute Diameter")
axs[2].set_xlabel("Parachute Diameter (m)")
axs[2].set_ylabel("Landing Success (1=Success, 0=Failure)")
axs[2].set_xticks(x)
axs[2].grid()

plt.tight_layout()
plt.savefig("parachute_drag_model_study.png")

