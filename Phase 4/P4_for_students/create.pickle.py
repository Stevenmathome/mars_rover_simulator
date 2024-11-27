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
 
 
 
planet = define_planet()
edl_system = define_edl_system()
mission_events = define_mission_events()
rover = define_rover()



df= pd.read_csv('rover_designs.csv')


df.drop(columns=['Unnamed: 0'], inplace=True) 

# print(df.columns)

min_index = df['total_time'].idxmin()
j = min_index


edl_system = redefine_edl_system(edl_system)
edl_system = define_batt_pack(edl_system,df.loc[j,'battery_type'], df.loc[j,'battery_number']) # battery type
edl_system = define_chassis(edl_system,df.loc[j,'chassis_type']) # chasis type
edl_system = define_motor(edl_system, df.loc[j,'motor_type']) # motor type


edl_system['parachute']['diameter'] = df.loc[j,'parachute_diameter']
edl_system['rover']['wheel_assembly']['wheel']['radius'] = df.loc[j,'wheel_radius']
edl_system['rover']['chassis']['mass'] =  df.loc[j,'chassis_mass']
edl_system['rover']['wheel_assembly']['speed_reducer']['diam_gear'] = df.loc[j,'gear_diameter']
edl_system['rocket']['initial_fuel_mass'] = df.loc[j,'fuel_mass']
edl_system['rocket']['fuel_mass'] = df.loc[j,'fuel_mass']

 
 
# tmax = 5000   
# experiment, end_event = experiment1()
# max_rover_velocity = -1  # this is during the landing phase
# min_strength=40000
# max_cost = 7.2e6
max_batt_energy_per_meter = edl_system['rover']['power_subsys']['battery']['capacity']/1000


 
edl_system['team_name'] = 'Pat the Bat'  # change this to something fun for your team (or just your team number)
edl_system['team_number'] = 17    # change this to your assigned team number (also change it below when saving your pickle file)



    # This will create a file that you can submit as your competition file.
with open('FA24_504team17.pickle', 'wb') as handle:
    pickle.dump(edl_system, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
# # Open the file in read-binary mode
with open('FA24_504team17.pickle', 'rb') as file:
    # Load the data from the file
    data = pickle.load(file)

print(data)