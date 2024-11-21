import pandas as pd
import numpy as np
import numpy as np
from subfunctions_Phase4 import *
from define_experiment import *
from scipy.optimize import minimize, differential_evolution
from scipy.optimize import Bounds
from scipy.optimize import NonlinearConstraint
import pickle
import sys
import pandas as pd
from csv import writer

df= pd.read_csv('rover_designs.csv')
df.drop(columns=['Unnamed: 0'], inplace=True) 

rows= 0 
drop = []
for j in range (df.shape[0]):
    
    planet = define_planet()
    edl_system = define_edl_system()
    mission_events = define_mission_events()
    rover = define_rover()




    def reset_initial_conditions(edl_system):
        edl_system['altitude'] = 11000    # [m] initial altitude
        edl_system['velocity'] = -587     # [m/s] initial velocity
        edl_system['parachute']['deployed'] = True   # our parachute is open
        edl_system['parachute']['ejected'] = False   # and still attached
        edl_system['rover']['on_ground'] = False # the rover has not yet landed
        
        
    battery = ['PbAcid-1', 'PbAcid-2', 'PbAcid-3', 'NiCD', 'NiMH','LiFePO4']
    motors  = ['base', 'base_he','torque', 'torque_he','speed', 'speed_he']



    edl_system = define_batt_pack(edl_system,df.iloc[j,3], df.iloc[j,4]) # battery type
    edl_system = define_chassis(edl_system,df.iloc[j,1]) # motor type
    edl_system = define_motor(edl_system, df.iloc[j,2]) # motor type


    
    tmax = 5000   
    experiment, end_event = experiment1()
    max_rover_velocity = -1  # this is during the landing phase
    min_strength=40000
    max_cost = 7.2e6
    max_batt_energy_per_meter = edl_system['rover']['power_subsys']['battery']['capacity']/1000



    reset_initial_conditions(edl_system)

    x= [df.iloc[j,5],df.iloc[j,6],df.iloc[j,7],df.iloc[j,8],df.iloc[j,9]]

    # check if we have a feasible solution 
    c = constraints_edl_system(x,edl_system,planet,mission_events,tmax,experiment,
                            end_event,min_strength,max_rover_velocity,max_cost,
                            max_batt_energy_per_meter)

    #c = [constraint_distance, constraint_strength, constraint_velocity, constraint_cost, constraint_battery]
    print(f"constraint_distance {c[0]:.1f}, constraint_strength {c[1]:.1f}, constraint_velocity {c[2]:.1f} constraint_cost {c[3]:.1f}, constraint_batter {c[4]:.1f}")
    print(np.max(c - np.zeros(len(c))))

    feasible = np.max(c - np.zeros(len(c))) <= 0
    if feasible:
        print('This is a feasible solution!')
        
    else:  # nonsense to let us know this did not work
        print('Solution not feasible')
        rows+=1
        drop.append(j)

df.drop(drop, inplace=True)
df.reset_index(drop=True, inplace=True)       

df.to_csv('rover_designs.csv')
print(f'{rows} rows were deleted')