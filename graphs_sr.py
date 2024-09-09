import matplotlib.pyplot as plt
import numpy as np


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
def get_gear_ratio(speed_reducer): # parker
  if (type(speed_reducer) != dict):
    raise Exception("The speed_reducr input must be a dictionary type")

  if speed_reducer['type'].lower() != 'reverted':
    raise Exception("The speed reducer you are using does not contain \'reverted'\ as its type")

  d2 = speed_reducer['diam_gear']
  d1 = speed_reducer['diam_pinion']
  Ng = (d2/d1)**2

  return Ng

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



motor = rover['wheel_assembly']['motor']
omega = np.array([1,2,3,4,5,6,7,8,9,10])
speed_reducer = rover['wheel_assembly']['speed_reducer']


Ng = get_gear_ratio(speed_reducer)
tau = tau_dcmotor(omega, motor)

if isinstance(tau, np.ndarray):
  tau_out = []
  omega_out = []
  power_out = []
  for i in range(len(tau)):
    tau_out.append(Ng * tau[i])
    omega_out.append(Ng * omega[i])
    power_out.append(tau_out[i] * omega_out[i])

  tau_out = np.array(tau_out)
  omega_out = np.array(omega_out)
  power_out = np.array(power_out)

if isinstance(tau, (int, float)):
  tau_out = tau * Ng
  omega_out = omega * Ng
  power_out = tau_out * omega_out

x = tau_out
y = omega_out

figs, axs = plt.subplots(nrows=3, ncols=1,figsize=(8, 12))
axs[0].scatter(x, y)
axs[0].set_title('motor shaft speed [rad/s] vs. motor shaft torque [Nm] ')
axs[0].set_xlabel("motor shaft torque [Nm]")
axs[0].set_ylabel("motor shaft speed [rad/s]")

x = tau_out
y = power_out

axs[1].scatter(x, y)
axs[1].set_title('motor power [W] vs. motor shaft torque [Nm]')
# axs[1].set_ylim(-100,0)
axs[1].set_xlabel("motor shaft torque [Nm]")
axs[1].set_ylabel("motor power [W]")

x = omega_out
y = power_out

axs[2].scatter(x, y)
axs[2].set_title('motor power [W] vs. motor shaft torque [Nm]')
# axs[1].set_ylim(-100,0)
axs[2].set_xlabel("motor shaft torque [Nm]")
axs[2].set_ylabel("motor power [W]")

plt.tight_layout()
plt.show()
