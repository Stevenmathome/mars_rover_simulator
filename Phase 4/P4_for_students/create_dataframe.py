import pandas as pd

df =  pd.DataFrame(columns=['chassis_type',
                            'battery_type',
                            'battery_number',    
                            'parachute_diameter',
                            'wheel_radius',
                            'chassis_mass',
                            'gear_diameter',
                            'fuel_mass',
                            'time_rover',
                            'total_time',
                            'edl_time',
                            'cost'
                            ])

df.to_csv('rover_designs1.csv')  

# df = pd.read_csv('rover_designs.csv',header = None, names = ['chassis_type',
#                             'battery_type',
#                             'battery_number',    
#                             'parachute_diameter',
#                             'wheel_radius',
#                             'chassis_mass',
#                             'gear_diameter',
#                             'fuel_mass',
#                             'time_rover',
#                             'total_time',
#                             'edl_time',
#                             'cost'
#                             ])

# df.to_csv('rover_designs.csv',index=False)


