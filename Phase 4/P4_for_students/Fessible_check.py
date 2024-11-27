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

from tqdm.auto import tqdm


df= pd.read_csv('rover_designs.csv')
# df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
rows= 0 
drop = []

for j in tqdm(range(df.shape[0]), desc="Processing rows", leave=True, dynamic_ncols=True, position=0, file=sys.stdout):
 # tqdm is a progress bar, it will show you how far along you are in the loop
    
    # if error is rasied delete the row
    try:
    
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



        edl_system = define_batt_pack(edl_system,df.loc[j,'battery_type'], df.loc[j,'battery_number']) # battery type
        edl_system = define_chassis(edl_system,df.loc[j,'chassis_type']) # chasis type
        edl_system = define_motor(edl_system, df.loc[j,'motor_type']) # motor type

        
        tmax = 5000   
        experiment, end_event = experiment1()
        max_rover_velocity = -1  # this is during the landing phase
        min_strength=40000
        max_cost = 7.2e6
        max_batt_energy_per_meter = edl_system['rover']['power_subsys']['battery']['capacity']/1000



        reset_initial_conditions(edl_system)

        x= [df.loc[j,'parachute_diameter'],df.loc[j,'wheel_radius'],df.loc[j,'chassis_mass'],df.loc[j,'gear_diameter'],df.loc[j,'fuel_mass']]

        # check if we have a feasible solution 
        c = constraints_edl_system(x,edl_system,planet,mission_events,tmax,experiment,
                                end_event,min_strength,max_rover_velocity,max_cost,
                                max_batt_energy_per_meter)

        #c = [constraint_distance, constraint_strength, constraint_velocity, constraint_cost, constraint_battery]
        # print(f"constraint_distance {c[0]:.1f}, constraint_strength {c[1]:.1f}, constraint_velocity {c[2]:.1f} constraint_cost {c[3]:.1f}, constraint_batter {c[4]:.1f}")
        # print(np.max(c - np.zeros(len(c))))

        feasible = np.max(c - np.zeros(len(c))) <= 0
        if feasible:
            # print('This is a feasible solution!')
            l=0
            
        else:   # nonsense to let us know this did not work
                 # print('Solution not feasible')
            rows+=1
            drop.append(j)

       
    
    except:
    #    drop.append(j)
    #    rows+=1
        print('nothing')

df.drop(drop, inplace=True)
df.reset_index(drop=True, inplace=True)

#find lowest 10 times and save it and get rid of the rest
df.sort_values(by='total_time', inplace=True)
df.reset_index(drop=True, inplace=True)
df = df.head(10)
       


df.to_csv('rover_designs.csv')
print(f'{rows} rows were deleted')