import matplotlib.pyplot as plt
import numpy as np
from subfunctions import tau_dcmotor
from subfunctions import get_gear_ratio


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
axs[0].plot(x, y)
axs[0].set_title('motor shaft speed [rad/s] vs. motor shaft torque [Nm] ')
axs[0].set_xlabel("motor shaft torque [Nm]")
axs[0].set_ylabel("motor shaft speed [rad/s]")

x = tau_out
y = power_out

axs[1].plot(x, y)
axs[1].set_title('motor power [W] vs. motor shaft torque [Nm]')
axs[1].set_xlabel("motor shaft torque [Nm]")
axs[1].set_ylabel("motor power [W]")

x = omega_out
y = power_out

axs[2].plot(x, y)
axs[2].set_title('motor power [W] vs. motor shaft torque [Nm]')
# axs[1].set_ylim(-100,0)
axs[2].set_xlabel("motor shaft torque [Nm]")
axs[2].set_ylabel("motor power [W]")

plt.tight_layout()
plt.show()