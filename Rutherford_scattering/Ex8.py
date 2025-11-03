# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 18:07:36 2025

@author: SS882604
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst

from my_package import Bellan_EMmotion as bemm

#%% INITIALIZATIONS

# ---------- LAB FRAME ----------

# LAB frame conditions
q1, q2 = -cst.e, -cst.e
m1, m2 = cst.m_e, cst.m_e

temp = np.loadtxt("Input_Ex8.txt")

init_Vel1, init_Vel2 = temp[0, :], temp[1, :]
init_Pos1, init_Pos2 = temp[2, :], temp[3, :]

# ---------- CM FRAME ----------

# NOTE: we are using the equivalent problem in the following computations
#
# Test particle: mass 'mu', charge 'q_test', velocity 'init_Vel', position 'init_Pos'
# Field (fixed) particle: infinite mass, charge 'q_field', null velocity, position in (0, 0)
#
# NOTE: r_relative = r_1 - r_2
# Thus: r_1 = r_cm + mu/m_1 * r_relative and r_2 = r_cm - mu/m_2 * r_relative

mu = m1 * m2 / (m1 + m2)
q_test = q1
q_field = q2

cm_Vel = (m1 * init_Vel1 + m2 * init_Vel2) / (m1 + m2)
cm_Pos = (m1 * init_Pos1 + m2 * init_Pos2) / (m1 + m2)
init_Vel = init_Vel1 - init_Vel2
init_Pos = init_Pos1 - init_Pos2

b_large = q_test * q_field / (4*cst.pi*cst.epsilon_0 * mu * np.dot(init_Vel, init_Vel))
b_large = np.abs(b_large)

print("\n-----------------------------------")
print("      Equivalent problem values     ")
print("-----------------------------------\n")

print(f"-----> Relative velocity: [{init_Vel[0]:8.2E}, {init_Vel[1]:8.2E}, {init_Vel[2]:8.2E}]")
print(f"-----> Relative position: [{init_Pos[0]:8.2E}, {init_Pos[1]:8.2E}, {init_Pos[2]:8.2E}]", end='\n\n')

# The 3D problem can be seen as casted in a 2D problem on the plane where:
#   - test particle is in 'init_Pos' and travels at velocity 'init_Vel'
#   - field particle is fixed in 0-0

# Given 'psi' as the angle between the relative velocity and position, I have that:
#   - init_Pos_mod * cos(psi): horizontal distance from the scattering center
#   - init_Pos_mod * sin(psi): impact parameter
#   - init_Vel_mod: the approaching speed (and velocity)
# All other entries are zero!

# I have empirically determined that: 
#   - 50 timesteps are enough to see the deflection
#   - the width of the time step is comparable to the collision time scale

init_Vel_mod = np.sqrt(np.dot(init_Vel, init_Vel))
init_Pos_mod = np.sqrt(np.dot(init_Pos, init_Pos))
factor = init_Vel_mod * init_Pos_mod

psi = np.arccos(np.dot(init_Vel, init_Pos)/factor)   # Approaching angle of the equivalent problem

tau = (b_large / init_Vel_mod) * 0.2   # Time scale of the collision. 1/5 is empirical
suggested_n_timesteps = int(110 + np.abs((init_Pos_mod * np.cos(psi) / init_Vel_mod) / tau))

print(f"-----> collision time-scale: {tau:.2E}")
print("-----> horizontal distance before collision: ", np.abs(init_Pos_mod * np.cos(psi)), end='\n\n')
print(f"IMPORTANT -----> large angle b = {b_large:.4E}")
print(f"IMPORTANT -----> impact parameter: {np.abs(init_Pos_mod * np.sin(psi)):.4E}")

print(f"\nSuggested distance before collision: {40 * b_large:.2E}")   # 40 is empirical
print("Suggested number of timesteps: ", suggested_n_timesteps)

B = np.array([0, 0, 0])
Omega = q_test * B / mu

#%% MAIN CODE

n_timesteps = 240
dt = tau   # For Rutherford scattering: 1e-7 with 150

# Copy without touching the initial values
Vel = np.copy(init_Vel)   # <--- relative velocity
Pos = np.copy(init_Pos)   # <--- relative position

all_Pos = np.array([Pos])   # Used for drawing the entire trajectory

for i in range(0, n_timesteps):
    
    E = bemm.compute_E(Pos, q_field)
    Sigma = q_test * E / mu
    
    New_Vel = bemm.compute_Vel(Vel, Omega, Sigma, dt)
    New_Pos = Pos + New_Vel * dt
    
    Vel = New_Vel
    Pos = New_Pos

    all_Pos = np.append(all_Pos, [Pos], axis=0)
    
# Positions in the CM frame
all_cmPos1 = all_Pos * mu / m1
all_cmPos2 = - all_Pos * mu / m2

# Positions in the LAB frame
time_axis = np.linspace(0, n_timesteps, n_timesteps + 1)*dt
all_cmPos = np.array([time_axis]*3)
all_cmPos = np.transpose(all_cmPos)
all_cmPos = all_cmPos * cm_Vel + cm_Pos

all_Pos1 = all_cmPos + all_cmPos1
all_Pos2 = all_cmPos + all_cmPos2

# ---------- PLOTS ----------

# bemm.plot_2DRutherford(all_Pos, all_cmPos1, all_cmPos2, all_Pos1, all_Pos2, cm_Vel)

bemm.plot_cmFrame(all_Pos)
bemm.plot_3DRutherford(all_cmPos1, all_cmPos2, "CM frame - Trajectories")
bemm.plot_3DRutherford(all_Pos1, all_Pos2, "LAB frame - Trajectories")

#%% BLA BLA

# Relative velocity and position, before and after rotation
init_Vel_mod = np.sqrt(np.dot(init_Vel, init_Vel))
init_Pos_mod = np.sqrt(np.dot(init_Pos, init_Pos))

alfa = - np.arctan2(init_Vel[1], init_Vel[0])
beta = - np.arcsin(init_Vel[1]/init_Vel_mod)

mat = bemm.give_RotMat('z', alfa)
rot_Vel = bemm.compute_MatVec(mat, init_Vel)
rot_Pos = bemm.compute_MatVec(mat, init_Pos)
mat = bemm.give_RotMat('y', beta)
rot_Vel = bemm.compute_MatVec(mat, rot_Vel)
rot_Pos = bemm.compute_MatVec(mat, rot_Pos)

versor = np.cross(rot_Vel, rot_Pos)/(init_Vel_mod*init_Pos_mod)
gamma = np.arcsin(versor[2])
mat = bemm.give_RotMat('x', gamma)
rot_Pos = bemm.compute_MatVec(mat, rot_Pos)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
bemm.plot_Vec(fig, ax, init_Vel/init_Vel_mod, 'blue', 'Vel - BEFORE')
bemm.plot_Vec(fig, ax, init_Pos/init_Pos_mod, 'green', 'Pos - BEFORE')
bemm.plot_Vec(fig, ax, rot_Vel/init_Vel_mod, 'aqua', 'Vel - AFTER')
bemm.plot_Vec(fig, ax, rot_Pos/init_Pos_mod, 'lime', 'Pos - AFTER')
ax.set_xlabel('X coord [m]')
ax.set_ylabel('Y coord [m]')
ax.set_zlabel('Z coord [m]')
ax.set_title("Relative vectors (normalized)")
ax.legend()
plt.show()

#%% 

# Thermal velocity: 1/2 * m * v**2 = 3/2 * k_B * T
T = 1   # In eV

v = np.sqrt(3 * cst.e * T / cst.m_e)   # v = 726392 m/s
l = np.cbrt(1/1e20)   # l = 2e-7 m
l
