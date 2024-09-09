# Parker
# plot three graphs 
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
tau = tau_dcmotor(omega, motor)


x=tau
y=omega

figs, axs = plt.subplots(nrows= 3 ,ncols=1,figsize=(8, 12))
axs[0].scatter(x,y)
axs[0].set_title('motor shaft speed [rad/s] vs. motor shaft torque [Nm] ')
axs[0].set_xlabel("motor shaft torque [Nm]")
axs[0].set_ylabel("motor shaft speed [rad/s]")


x= tau
y= []
for i in range(len(omega)):
  y.append(omega[i]*tau[i])
y = np.array(y)

axs[1].scatter(x,y)
axs[1].set_title('motor power [W] vs. motor shaft torque [Nm]')
axs[1].set_xlabel("motor shaft torque [Nm]")
axs[1].set_ylabel("motor power [W]")



x= omega 
y= []
for i in range(len(omega)):
  y.append(omega[i]*tau[i])
y = np.array(y)

axs[2].scatter(x,y)
axs[2].set_title('motor power [W] vs. motor shaft speed [rad/s]')

axs[2].set_xlabel("motor shaft speed [rad/s] ")
axs[2].set_ylabel("motor power [W]")

plt.tight_layout()
plt.show()
