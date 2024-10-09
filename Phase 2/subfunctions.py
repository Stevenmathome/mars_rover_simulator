import numpy as np
import math
from scipy.interpolate import interp1d
from Utils.define_experiment import experiment1
import scipy.integrate as spi
# Constants

# Dictionaries
# this has upaded rover[‘wheel_assembly’][‘motor’]
rover = {
  'wheel_assembly': {
      'wheel': {'radius': 0.30, 'mass': 1.0
                },
      
      'motor': {'torque_stall': 170, 'torque_noload': 0, 'speed_noload': 3.80, 'mass': 5.0, 
                'effcy_tau' : np.array([0, 10, 20, 40, 70, 165]) ,
     'effcy' : np.array([0, .55,.75,.71,.50,.05])
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


experiment , end_event = experiment1()
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
    raise Exception("The motor input must be a dictionary type")
  elif isinstance(omega, np.ndarray):
    if omega.ndim != 1:
        raise Exception("First input omega must be a scalar or vector numpy array. Matrices are not allowed")

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
      if omega.ndim != 1:
        raise Exception("The omega input must be a 1D numpy array")
      if  omega[i] > motor["speed_noload"]:
        tau.append(float(motor["torque_noload"]))
      elif omega[i] < 0 :
        tau.append(float(motor['torque_stall']))
      elif  0 <= omega[i] <= motor["speed_noload"]:
        tau.append(float (motor['torque_stall'] - ((motor['torque_stall'] - motor["torque_noload"])/motor["speed_noload"])*omega[i]))

    tau = np.array(tau)

    return tau


def F_drive(omega,rover): 
  
  if (type(rover) != dict):
    raise Exception("The second input rover must be a dictionary.")
  if rover is None:
        raise ValueError("The rover input cannot be None")
    # Check if omega is None
  if omega is None:
        raise ValueError("The omega input cannot be None")
  #elif isinstance(omega, (int, float)):
  #  if np.isnan(omega) or np.isinf(omega):
  #    raise Exception("The omega input cannot contain NaN or Inf values")
  if isinstance(omega, list):
    raise Exception("The first input 'omega' must be a scalar or numpy array.")
  if isinstance(omega, np.ndarray):
        # Check if omega contains None values
        if np.any(omega == None):
            raise ValueError("The omega input cannot contain None values")
        # Check if omega contains NaN or Inf
        if np.any(np.isnan(omega)):
            raise ValueError("The omega input cannot contain NaN values")
        if np.any(np.isinf(omega)):
            raise ValueError("The omega input cannot contain Inf values")
        if omega.ndim != 1:
            raise ValueError

          
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
        raise Exception("Terrain angle must be a scalar or a 1D numpy array")

    if isinstance(terrain_angle, np.ndarray) and not np.all((-75 <terrain_angle) & (terrain_angle < 75)):
        raise Exception("Terrain angle must be between -75 and 75 degrees")

    if not isinstance(rover, dict):
        raise Exception("Must be dict")

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
    radius = rover['wheel_assembly']['wheel']['radius']
    # Convert angle to radians
    if isinstance(terrain_angle,(int,float)):
      terrain_rad = np.radians(terrain_angle)
      # Compute normal force
      N_force = (mass * g * np.cos(terrain_rad)) / 6
      # Compute rolling resistance force
      wheel_omega = omega / gear_ratio
      frr = (-Crr * N_force * math.erf(40*radius*wheel_omega)) * 6
      return frr
      
    else:
      frr=[]
      for i in range(len(terrain_angle)):
        terrain_rad = np.radians(terrain_angle[i])
        # Compute normal force
        N_force = (mass * g * np.cos(terrain_rad)) / 6
        # Compute rolling resistance force
        wheel_omega = omega[i] / gear_ratio
        frr.append((-Crr * N_force * math.erf(40*radius*wheel_omega)) * 6)
      frr = np.array(frr)
      return frr

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

def motorW(v, rover): 
    if not isinstance(v, (int,float,np.ndarray)):
        raise Exception ("v must be and scalar or 1d array")
    if np.ndim(v)!= 1:
        raise Exception('v must be a 1d array')
    if not isinstance (rover, (dict)):
        raise Exception("Rover must be a dictionary")
    
    ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

    
    w= (v / rover['wheel_assembly']['wheel']['radius'])/ng
    
    return w
    #inputs should be 1D array/ scalar and a dictionary\
        
def rover_dynamics(t, y, rover, planet, experiment):
    if not isinstance(t , (int , float)):
        raise Exception(" t must be a scalar")
    if not isinstance (y, (np.ndarray)):
        raise Exception(" y must be a numpy array")
    if np.ndim(y)!=1:
        raise Exception("y must be a 1D array with two elements")
    if len(y)!=2:
        raise Exception("y must have two elements")
    if not isinstance(rover, (dict)):
        raise Exception("rover must be a dictionary")
    if not isinstance(planet, (dict)):
        raise Exception("planet must be a dictionary")
    if not isinstance(experiment, (dict)):
        raise Exception("experiment must be a dictionary")
    
    alpha_dist = experiment['alpha_dist']   
    alpha_deg = experiment['alpha_deg']
    
    alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate')
    # fit the cubic spline
    terrian_angle = alpha_fun(y[1])
    w = motorW(y[0], rover)
    
    F_net1 = F_net(w, terrian_angle, rover, planet, experiment['Crr'])
    
    acceleration = F_net1 / get_mass(rover)
    
    velocity = y[0] + acceleration * t
    
    return np.array([acceleration, velocity])


def mechpower(v,rover):
    if not isinstance(v, (int,float,np.ndarray)):
        raise Exception ("v must be and scalar or 1d array")
    if np.ndim(v)!= 1:
        raise Exception('v must be a 1d array')
    if not isinstance (rover, (dict)):
        raise Exception("Rover must be a dictionary")
    
    w = motorW(v, rover)
    
    tau = tau_dcmotor(w, rover['wheel_assembly']['motor'])
    
    p = tau * w
    
    return p
    #inputs should be 1D array/ scalar and a dictionary

def battenergy(t,v,rover):
    if not isinstance(t,np.ndarray):
        raise Exception("t must be a numpy array")
    if not isinstance(v,np.ndarray):
        raise Exception("v must be a numpy array")
    if len(t)!=len(v):
        raise Exception("t and v must have the same length")
    if not isinstance(rover,dict):
        raise Exception("rover must be a dictionary")
    if np.ndim(t)!=1:
        raise Exception("t must be a 1D array")
    if np.ndim(v)!=1:
        raise Exception("v must be a 1D array")
    
    tau = tau_dcmotor(v,rover['wheel_assembly']['motor'])
    p = mechpower(v,rover)
    effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
    effcy = rover['wheel_assembly']['motor']['effcy']
    
    effcy_fun = interp1d(effcy_tau, effcy, kind = 'cubic', fill_value='extrapolate')
    # x = np.linspace(tau.min, tau.max, 100)
    
    # print(p)
    # print(effcy_fun(tau))
    pbatt = p / effcy_fun(tau)
  
    E = spi.trapezoid(pbatt, t) * 6
    
    return E
    
def simulate_rover(rover,planet,experiment,end_event):
    if not isinstance(rover,dict):
        raise Exception("rover must be a dictionary")
    if not isinstance(planet,dict):
        raise Exception("planet must be a dictionary")
    if not isinstance(experiment,dict):
        raise Exception("experiment must be a dictionary")
    if not callable(end_event):
        raise Exception("end_event must be a function")