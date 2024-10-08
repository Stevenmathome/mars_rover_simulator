# Parker
# plot three graphs 
import matplotlib.pyplot as plt
import numpy as np
from subfunctions import tau_dcmotor

  
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
omega = np.linspace(0,3.8,100)
tau = tau_dcmotor(omega, motor)

x=tau
y=omega

figs, axs = plt.subplots(nrows= 3 ,ncols=1,figsize=(8, 15))
axs[0].plot(x,y)
axs[0].set_title('motor shaft speed [rad/s] vs. motor shaft torque [Nm] ')
axs[0].set_xlabel("motor shaft torque [Nm]")
axs[0].set_ylabel("motor shaft speed [rad/s]")


x= tau
y= []
for i in range(len(omega)):
  y.append(omega[i]*tau[i])
y = np.array(y)

axs[1].plot(x,y)
axs[1].set_title('motor power [W] vs. motor shaft torque [Nm]')
axs[1].set_xlabel("motor shaft torque [Nm]")
axs[1].set_ylabel("motor power [W]")



x= omega 
y= []
for i in range(len(omega)):
  y.append(omega[i]*tau[i])
y = np.array(y)

axs[2].plot(x,y)
axs[2].set_title('motor power [W] vs. motor shaft speed [rad/s]')
axs[2].set_xlabel("motor shaft speed [rad/s] ")
axs[2].set_ylabel("motor power [W]")

plt.tight_layout()
plt.show()