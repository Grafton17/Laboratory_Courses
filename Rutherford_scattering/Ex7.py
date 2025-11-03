# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
import matplotlib as mpl

from my_package import Bellan_EMmotion as bemm
    
#%% INITIALIZATIONS

# Strange behaviour in the circular motion when only B is applied, the particle
# seems skipping one half of the circle in its trajectory...

# Boundary conditions
q = + cst.e
m = cst.m_e
B = np.array([0, 0, 0])   # In Teslas
E = np.array([1, 0, 0])   # In Volts/meter

Omega = q*B/m
Sigma = q*E/m

# Initial conditions
init_Vel = np.array([1000, 0, 0])
init_Pos = np.array([-0.02, 0.001, 0])

# The time step is determined by the boundary conditions:
#   - accelerating motion along the E line
#   - cyclotron motion around the B line

# ELECTRIC MOTION:
#   a = qE/m

# CYCLOTRON MOTION:
#   r = mv/qB sets the length scale
#   w = qB/m or t ~ m/qB sets the time scale
#   With an E field, we have: v ~ a*t ~ E/B

B_module = np.sqrt(np.dot(B, B))
E_module = np.sqrt(np.dot(E, E))
w = q*B_module / m
acc = q*E_module / m

print("\n-------------------------------------")
print("      Set the suitable timestep     ")
print("-------------------------------------\n")

if (B_module != 0):
    print(f"-----> Cyclotron motion detected:", end="\t")
    print(f"Frequency: {w:.2e} - Time scale: {2*cst.pi/w:.2e}")
else:
    print(f"-----> No cylotron motion detected:", end="\t")
    print(f"Acceleration: {acc:.2e} - Time per 1 meter: {np.sqrt(2/acc):.2e}")

#%% MAIN CODE

dt = 3 * 1e-7   # For cyclotron motion:at least 0.01 * 1/cyclotron_frequency
n_timesteps = 150

# Copy without touching the initial values
Pos = np.copy(init_Pos)
Vel = np.copy(init_Vel)

all_Pos = np.array([Pos])   # Used for drawing the entire trajectory

for i in range(0, n_timesteps):
    
    E = bemm.compute_E(Pos, cst.e)
    Sigma = cst.e * E / cst.m_e
    
    New_Vel = bemm.compute_Vel(Vel, Omega, Sigma, dt)
    New_Pos = Pos + New_Vel * dt
    
    Vel = New_Vel
    Pos = New_Pos

    all_Pos = np.append(all_Pos, [Pos], axis=0)
#%% PLOT

fig, [axP, axT] = plt.subplots(2, 1)

axP.plot(all_Pos[:, 0], all_Pos[:, 1], '.')
axP.axis('equal')
axP.set_xlabel("x position [m]")
axP.set_ylabel("y position [m]") 
axP.set_title("Trajectory")
axP.grid()

time_axis = np.linspace(0, n_timesteps, n_timesteps + 1)*dt

time_law = bemm.nothing(time_axis)
# time_law = bemm.MRUA(time_axis, acc)
# time_law = bemm.HM(time_axis, m * init_Vel[0] / q * B[2], w)
# time_law = bemm.ExB(time_axis, E_module/B_module, w)

axT.plot(time_axis, all_Pos[:, 1], label='Simulated')
axT.plot(time_axis, time_law, label='Expected')
axT.set_xlabel("Time [s]")
axT.set_ylabel("y position [m]")
axT.set_title("Time law")
axT.grid()
axT.legend()

plt.tight_layout()
plt.show()
#%% BRUTTA

cst.e
