
"""
Created on Mon Nov 15 22:58:03 2021

@author: Marvin Engineering Design Team
"""

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

# the following calls instantiate the needed structs and also make some of
for i in range(5):
    # our design selections (battery type, etc.)
    planet = define_planet()
    edl_system = define_edl_system()
    mission_events = define_mission_events()
    rover = define_rover()
    df= pd.read_csv('rover_designs.csv')

    # all possible things to define in the edl_system

    battery = ['PbAcid-1', 'PbAcid-2', 'PbAcid-3', 'NiCD', 'NiMH','LiFePO4']
    motors  = ['base', 'base_he','torque', 'torque_he','speed', 'speed_he']

    #constraints
    edl_system = define_chassis(edl_system,'magnesium') # motor type
    edl_system = define_motor(edl_system, motors[2]) # motor type
    edl_system = define_batt_pack(edl_system,battery[2], 10)

    tmax = 5000

    # Overrides what might be in the loaded data to establish our desired
    # initial conditions
    def reset_initial_conditions(edl_system):
        edl_system['altitude'] = 11000    # [m] initial altitude
        edl_system['velocity'] = -587     # [m/s] initial velocity
        edl_system['parachute']['deployed'] = True   # our parachute is open
        edl_system['parachute']['ejected'] = False   # and still attached
        edl_system['rover']['on_ground'] = False # the rover has not yet landed
        
        
        return edl_system

    reset_initial_conditions(edl_system)

    experiment, end_event = experiment1()

    # constraints
    max_rover_velocity = -1  # this is during the landing phase
    min_strength=40000
    max_cost = 7.2e6
    max_batt_energy_per_meter = edl_system['rover']['power_subsys']['battery']['capacity']/1000


    # ******************************
    # DEFINING THE OPTIMIZATION PROBLEM
    # ****
    # Design vector elements (in order):
    #   - parachute diameter [m]
    #   - wheel radius [m]
    #   - chassis mass [kg]
    #   - speed reducer gear diameter (d2) [m]
    #   - rocket fuel mass [kg]
    #

    # search bounds
    #x_lb = np.array([14, 0.2, 250, 0.05, 100])
    #x_ub = np.array([19, 0.7, 800, 0.12, 290])
    bounds = Bounds([17, 0.3, 250, 0.05, 100], [19, 0.7, 800, 0.12, 290])

    # initial guess
    # x0 = np.array([np.random.uniform(17,19), np.random.uniform(0.3,0.7), np.random.uniform(250,800), np.random.uniform(0.05,0.12), np.random.uniform(100,290)])
    # x0 = np.array([19, .7, 550.0, 0.09, 250.0]) 

    # lambda for the objective function
    obj_f = lambda x: obj_fun_time(x,
                                edl_system,
                                planet,
                                mission_events,
                                tmax,
                                experiment,end_event
                                )

    # lambda for the constraint functions
    #   ineq_cons is for SLSQP
    #   nonlinear_constraint is for trust-constr
    cons_f = lambda x: constraints_edl_system(x,edl_system,planet,mission_events,
                                            tmax,experiment,end_event,min_strength,
                                            max_rover_velocity,max_cost,max_batt_energy_per_meter)
    nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 0)  # for trust-constr
    ineq_cons = {'type' : 'ineq',
                'fun' : lambda x: -1*constraints_edl_system(x,edl_system,planet,
                                                            mission_events,tmax,experiment,
                                                            end_event,min_strength,max_rover_velocity,
                                                            max_cost,max_batt_energy_per_meter)}
    Nfeval = 1
    def callbackF(Xi):  # this is for SLSQP reporting during optimization
        global Nfeval
        if Nfeval == 1:
            print('Iter    Par_Diam     Wheel Rad    Chassis Mass  Gear Diam  Fuel Mass        fval')
            
        print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}  {5: 3.6f} \
            {6: 3.6f}'.format(Nfeval, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], obj_f(Xi)))
        Nfeval += 1



    # The optimizer options below are
    # 'trust-constr'
    # 'SLSQP'
    # 'differential_evolution'
    # 'COBYLA'
    # You should fully comment out all but the one you wish to use

    ###############################################################################
    #call the trust-constr optimizer --------------------------------------------#
    # options = {'maxiter': 5, 
    #             # 'initial_constr_penalty' : 5.0,
    #             # 'initial_barrier_parameter' : 1.0,
    #             'verbose' : 3,
    #             'disp' : True}
    # res = minimize(obj_f, x0, method='trust-constr', constraints=nonlinear_constraint, 
    #                 options=options, bounds=bounds)
    # end call to the trust-constr optimizer -------------------------------------#
    ###############################################################################

    ###############################################################################
    # call the SLSQP optimizer ---------------------------------------------------#
    #testing different mass chasis

    edl_system['rover']['chassis']['mass'] = 250
    mass= get_mass_edl(edl_system)
    strength = 250* edl_system['rover']['chassis']['mass']
    print(f"mass {mass:.1f}, strength {strength:.1f}")

    # min_mass_steel = 40000
    # min_mass_magnesium
    # min_mass_carbon_fiber


    res_save = []
    for i in range(2):
        x0 = np.array(
            [np.random.uniform(14,19), #parachute diameter
            np.random.uniform(0.2,0.7), #wheel radius
            # np.random.uniform(250,800), #chassis mass
            250, #chassis mass
            np.random.uniform(0.05,0.12), #speed reducer gear diameter
            np.random.uniform(100,290)] #rocket fuel mass
                    )
        
        reset_initial_conditions(edl_system)
        options = {'maxiter': 10,
                    'disp' : True}
        res = minimize(obj_f, x0, method='SLSQP', constraints=ineq_cons, bounds=bounds, 
                        options=options, callback=callbackF)
        
        # check if we have a feasible solution 
        c = constraints_edl_system(res.x,edl_system,planet,mission_events,tmax,experiment,
                                end_event,min_strength,max_rover_velocity,max_cost,
                                max_batt_energy_per_meter)
        
        #c = [constraint_distance, constraint_strength, constraint_velocity, constraint_cost, constraint_battery]
        print(f"constraint_distance {c[0]:.1f}, constraint_strength {c[1]:.1f}, constraint_velocity {c[2]:.1f} constraint_cost {c[3]:.1f}, constraint_batter {c[4]:.1f}")
        print(np.max(c - np.zeros(len(c))))
        
        feasible = np.max(c - np.zeros(len(c))) <= 0
        if feasible:
            xbest = res.x
            fbest = res.fun
            res_save.append(res.x)
        else:  # nonsense to let us know this did not work
            xbest = [99999, 99999, 99999, 99999, 99999]
            fval = [99999]
            print('Solution not feasible, restarting tests...')
        Nfeval = 1
        
    if len(res_save) == 0:
        raise Exception('Solution not feasible, exiting code...')
    # end call to the SLSQP optimizer --------------------------------------------#
    ###############################################################################

    ###############################################################################
    # call the differential evolution optimizer ----------------------------------#
    # popsize=2 # define the population size
    # maxiter=20 # define the maximum number of iterations
    # res = differential_evolution(obj_f, bounds=bounds, constraints=nonlinear_constraint, popsize=popsize, maxiter=maxiter, disp=True, polish = False) 
    # end call the differential evolution optimizer ------------------------------#
    ###############################################################################

    ###############################################################################
    # call the COBYLA optimizer --------------------------------------------------#
    # cobyla_bounds = [[14, 19], [0.2, 0.7], [250, 800], [0.05, 0.12], [100, 290]]
    # #construct the bounds in the form of constraints
    # cons_cobyla = []
    # for factor in range(len(cobyla_bounds)):
    #     lower, upper = cobyla_bounds[factor]
    #     l = {'type': 'ineq',
    #           'fun': lambda x, lb=lower, i=factor: x[i] - lb}
    #     u = {'type': 'ineq',
    #           'fun': lambda x, ub=upper, i=factor: ub - x[i]}
    #     cons_cobyla.append(l)
    #     cons_cobyla.append(u)
    #     cons_cobyla.append(ineq_cons)  # the rest of the constraints
    # options = {'maxiter': 50, 
    #             'disp' : True}
    # res = minimize(obj_f, x0, method='COBYLA', constraints=cons_cobyla, options=options)
    # end call to the COBYLA optimizer -------------------------------------------#
    ###############################################################################


    # # check if we have a feasible solution 
    # c = constraints_edl_system(res.x,edl_system,planet,mission_events,tmax,experiment,
    #                            end_event,min_strength,max_rover_velocity,max_cost,
    #                            max_batt_energy_per_meter)

    # feasible = np.max(c - np.zeros(len(c))) <= 0
    # if feasible:
    #     xbest = res.x
    #     fbest = res.fun
    # else:  # nonsense to let us know this did not work
    #     xbest = [99999, 99999, 99999, 99999, 99999]
    #     fval = [99999]
    #     raise Exception('Solution not feasible, exiting code...')
    #     sys.exit()

    # What about the design variable bounds?

    # The following will rerun your best design and present useful information
    # about the performance of the design
    # This will be helpful if you choose to create a loop around your optimizers and their initializations
    # to try different starting points for the optimization.
    total_cost = []
    for i in range(len(res_save)):
        xbest = res_save[i]
        edl_system = redefine_edl_system(edl_system)

        edl_system['parachute']['diameter'] = xbest[0]
        edl_system['rover']['wheel_assembly']['wheel']['radius'] = xbest[1]
        edl_system['rover']['chassis']['mass'] = xbest[2]
        edl_system['rover']['wheel_assembly']['speed_reducer']['diam_gear'] = xbest[3]
        edl_system['rocket']['initial_fuel_mass'] = xbest[4]
        edl_system['rocket']['fuel_mass'] = xbest[4]
        total_cost.append(get_cost_edl(edl_system))

    ind = total_cost.index(min(total_cost)) 
    xbest = res_save[ind]  

    # *****************************************************************************
    # These lines save your design for submission for the rover competition.
    # You will want to change them to match your team information.

    edl_system['team_name'] = 'FunTeamName'  # change this to something fun for your team (or just your team number)
    edl_system['team_number'] = 4    # change this to your assigned team number (also change it below when saving your pickle file)

    # This will create a file that you can submit as your competition file.
    with open('FA24_501team99.pickle', 'wb') as handle:
        pickle.dump(edl_system, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # *****************************************************************************

    #del edl_system
    #with open('challenge_design_team9999.pickle', 'rb') as handle:
    #    edl_system = pickle.load(handle)

    time_edl_run,_,edl_system = simulate_edl(edl_system,planet,mission_events,tmax,True)
    time_edl = time_edl_run[-1]

    edl_system['rover'] = simulate_rover(edl_system['rover'],planet,experiment,end_event)
    time_rover = edl_system['rover']['telemetry']['completion_time']

    total_time = time_edl + time_rover
    
    edl_system_total_cost=get_cost_edl(edl_system)

    print('----------------------------------------')
    print('----------------------------------------')
    print('Optimized parachute diameter   = {:.6f} [m]'.format(xbest[0]))
    print('Optimized rocket fuel mass     = {:.6f} [kg]'.format(xbest[4]))
    print('Time to complete EDL mission   = {:.6f} [s]'.format(time_edl))
    print('Rover velocity at landing      = {:.6f} [m/s]'.format(edl_system['rover_touchdown_speed']))
    print('Optimized wheel radius         = {:.6f} [m]'.format(xbest[1])) 
    print('Optimized d2                   = {:.6f} [m]'.format(xbest[3])) 
    print('Optimized chassis mass         = {:.6f} [kg]'.format(xbest[2]))
    print('Time to complete rover mission = {:.6f} [s]'.format(time_rover))
    print('Time to complete mission       = {:.6f} [s]'.format(total_time))
    print('Average velocity               = {:.6f} [m/s]'.format(edl_system['rover']['telemetry']['average_velocity']))
    print('Distance traveled              = {:.6f} [m]'.format(edl_system['rover']['telemetry']['distance_traveled']))
    print('Battery energy per meter       = {:.6f} [J/m]'.format(edl_system['rover']['telemetry']['energy_per_distance']))
    print('Total cost                     = {:.6f} [$]'.format(edl_system_total_cost))
    print('----------------------------------------')
    print('----------------------------------------')

    new_row = {
        'chassis_type': edl_system['rover']['chassis']['type'],
        'battery_type': edl_system['rover']['power_subsys']['battery']['battery_type'],
        'battery_number': edl_system['rover']['power_subsys']['battery']['num_modules'],
        'parachute_diameter': xbest[0],
        'wheel_radius': xbest[1],
        'chassis_mass': xbest[2],
        'gear_diameter': xbest[3],
        'fuel_mass': xbest[4],
        'total_time': total_time,
        "edl_time": time_edl,
        'time_rover': time_rover,
        'cost': total_cost[ind],
        'motor_type': edl_system['rover']['wheel_assembly']['motor']['type']
    }


    # Convert the new row to a DataFrame and concatenate it with the existing one
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    print(df)

    # Append the new DataFrame row to the CSV file
    df.to_csv('rover_designs.csv', mode='a', index=False, header=False)