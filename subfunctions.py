import numpy as np
# Constants

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

# Functions
def get_mass(rover): # parker
  if (type(rover) != dict):
    raise Exception("The rover input must be a dictionary type")

  mass = 0

  mass += 6 * float(rover['wheel_assembly']['wheel']['mass'])
  mass += 6 * float(rover['wheel_assembly']['speed_reducer']['mass'])
  mass += 6 * float(rover['wheel_assembly']['motor']['mass'])

  mass += rover['chassis']['mass']

  mass += rover["science_payload"]['mass']

  mass += rover["power_subsys"]['mass']

  return mass



def get_gear_ratio(speed_reducer): # parker
  if (type(speed_reducer) != dict):
    raise Exception("The speed_reducr input must be a dictionary type")

  if speed_reducer['type'].lower() != 'reverted':
    raise Exception("The speed reducer you are using does not contain reverted as its type")

  d2 = speed_reducer['diam_gear']
  d1 = speed_reducer['diam_pinion']
  Ng = (d2/d1)**2

  return Ng


def tau_dcmotor(omega, motor): # parker

  if (type(motor) != dict):
    raise Exception("The omega input must be a dictionary type")
  if not isinstance(omega, (int, float, np.ndarray)):
      raise Exception("The omega input must be of type int, float, or np.ndarray")

  if isinstance(omega, (int, float)):
    tau = 0
    if  omega > motor["speed_noload"]:
      tau  = motor["torque_noload"]
    if  omega < 0 :
      tau = motor['torque_stall']
    if  0 <= omega <= motor["speed_noload"]:
      tau = motor['torque_stall'] - ((motor['torque_stall'] - motor["torque_noload"])/motor["speed_noload"])*omega

    return tau

  if isinstance(omega, np.ndarray):
    tau =  []
    for i in range (len(omega)):
      if  omega[i] > motor["speed_noload"]:
        tau.append(float(motor["torque_noload"]))
      if  omega[i] < 0 :
        tau.append(float(motor['torque_stall']))
      if  0 <= omega[i] <= 3.8:
        tau.append(float (motor['torque_stall'] - ((motor['torque_stall'] - motor["torque_noload"])/motor["speed_noload"])*omega[i]))

    tau = np.array(tau)

    return tau


def F_drive(omega,rover): 
    if not (np.isscalar(omega) or (isinstance(omega,np.ndarray) and omega.ndim in [0,1])): 
        raise Exception("Speical Message")

    if not isinstance(rover,dict):
        raise Exception ("Special ")

    motor = rover['wheel_assembly']['motor']
    speed_reducer = rover['wheel_assembly']['speed_reducer']
    wheel = rover['wheel_assembly']['wheel']

    tau_motor = tau_dcmotor(omega, motor)
    N_gear = get_gear_ratio(speed_reducer)
    
    tau_wheel = N_gear * tau_motor
    F_wheel = tau_wheel / wheel['radius']
    
    return 6 * F_wheel  # Multiply by 6 for all wheels  

def F_gravity(terrain_angle, rover, planet): 
    if not (np.isscalar(terrain_angle) or (isinstance(terrain_angle, np.ndarray) and terrain_angle.ndim in [0, 1])):
        raise Exception("Special Message")

    if isinstance(terrain_angle, np.ndarray) and not np.all((-75 <= terrain_angle) & (terrain_angle <= 75)):
        raise Exception("Special Message ")

    if not isinstance(rover, dict):
        raise Exception("Must be dict ")

    if not isinstance(planet, dict):
        raise Exception("must be dict ")

    # Extract parameters
    m = get_mass(rover)  
    g = planet['g']  

    # Convert terrain angle from degrees to radians
    theta_rad = np.radians(terrain_angle)

    # Compute the gravitational force component in the direction of           rover translation
    Fgt = -m * g * np.sin(theta_rad)

    return Fgt


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    if not (np.isscalar(omega) and np.isscalar(terrain_angle)) and not (isinstance(omega, np.ndarray) and isinstance(terrain_angle, np.ndarray) and omega.shape == terrain_angle.shape):
        raise Exception("The first two inputs must be either scalars or numpy arrays of the same shape.")

    if np.any(np.abs(terrain_angle) > 75):
        raise Exception("All elements of terrain_angle must be between -75 degrees and +75 degrees.")

    if not isinstance(rover, dict):
        raise Exception("The rover input must be a dictionary type")

    if not isinstance(planet, dict):
        raise Exception("The planet input must be a dictionary type")

    if not isinstance(Crr, (int, float)) or Crr <= 0:
        raise Exception("The coefficient of rolling resistance (Crr) must be a positive scalar")

    mass = get_mass(rover)
    g = planet['g']
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

    # Convert angle to radians
    terrain_rad = np.radians(terrain_angle)

    # Compute normal force
    N_force = (mass * g * np.cos(terrain_rad)) / 6

    # Compute rolling resistance force
    wheel_omega = omega * gear_ratio
    return -Crr * N_force * np.sign(wheel_omega) * 6

def F_net(omega, terrain_angle, rover, planet, Crr): #steve
  # Input validation
  if not (np.isscalar(omega) and np.isscalar(terrain_angle)) and not (isinstance(omega, np.ndarray) and isinstance(terrain_angle, np.ndarray) and omega.shape == terrain_angle.shape):
    raise Exception("The first two inputs must be either scalars or numpy arrays of the same shape.")


  if np.any(np.abs(terrain_angle) > 75):
    raise Exception("All elements of terrain_angle must be between -75 degrees and +75 degrees.")

  if (type(rover) != dict):
    raise Exception("The rover input must be a dictionary type")

  if (type(planet) != dict):
    raise Exception("The planet input must be a dictionary type")

  if not isinstance(Crr, (int, float)) or Crr <= 0:
    raise Exception("The coefficient of rolling resistance (Crr) must be a positive scalar")


  F_d = F_drive(omega, rover)
  F_g = F_gravity(terrain_angle, rover, planet)
  F_r = F_rolling(omega, terrain_angle, rover, planet, Crr)
    
  return F_d + F_g + F_r  