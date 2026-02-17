# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 13:11:03 2025

@author: SS882604
"""

folder_path = r"Y:\My Drive\MIA_UNI\Magistrale\CMS\Working_folder\L20-22_03Dicembre_exercise"

import numpy as np
import matplotlib.pyplot as plt
import my_package.my_module as md
import matplotlib as mpl

#%% FAST CHECK

output_folder = folder_path + "/" + "point5_1300K_2fs_300ps"   # point4_1700K_2fs_300ps

data = np.loadtxt(output_folder + "/Output.txt")
time, evolution_Etot, evolution_T = data[:,1], data[:, 2], data[:, 5]

md.compute_RelFluct(evolution_Etot, time, verbose=True, thermal_time = 10 * 1e12)
md.compute_RelFluct(evolution_T, time, verbose=True, thermal_time = 10 * 1e12)

# https://stackoverflow.com/questions/17788685/python-saving-multiple-figures-into-one-pdf-file

fig, axE, axT = md.plot_timeEvolution(time, evolution_Etot, evolution_T)



# When you need to compare two different simulations

# output_folder = folder_path + "/" + "noPOLY_500K_1fs_5000iter"

# data = np.loadtxt(output_folder + "/Output.txt")
# time, evolution_Etot, evolution_T = data[:,1], data[:, 2], data[:, 5]

# md.compute_RelFluct(evolution_Etot, time, verbose=True)

# axE.plot(time, evolution_Etot)
# axT.plot(time, evolution_T)




del output_folder
del data
del time
del evolution_Etot
del evolution_T
del fig, axE, axT

#%% EXERCISE 1.1 and 1.2 ----- ANALYZE DATA

sim_names = [ 
    "POLY_2000K_1fs_5000iter",
    "POLY_2000K_2fs_5000iter",
    "POLY_2000K_3fs_5000iter",
    "POLY_2000K_4fs_5000iter",
    "POLY_2000K_5fs_5000iter",
    "POLY_2000K_6fs_5000iter"  
    ]

sim_dt = [1, 2, 3, 4, 5, 6]

with open(folder_path + "/1_results_POLY_2000K.txt", "w") as f_results:

    for sim_name, dt in zip(sim_names, sim_dt):
        output_folder = folder_path + "/" + sim_name
        
        try:
            data = np.loadtxt(output_folder + "/Output.txt")
        except FileNotFoundError:
            del sim_names, sim_dt, sim_name, dt
            del output_folder, f_results
            del data
            del time, evolution_Etot, evolution_T
            del E_mean, E_fluct
            del T_mean, T_fluct
            del dummy
            raise
    
        time, evolution_Etot, evolution_T = data[:,1], data[:, 2], data[:, 5]

        E_mean, dummy, E_fluct = md.compute_RelFluct(evolution_Etot, time, verbose=False)
        T_mean, dummy, T_fluct = md.compute_RelFluct(evolution_T, time, verbose=False)
    
        print(sim_name)
        f_results.write(f"{str(dt).zfill(2)}\t{E_mean}\t{E_fluct}\t{T_mean}\t{T_fluct}\n")
        
del sim_names, sim_dt, sim_name, dt
del output_folder, f_results
del data
del time, evolution_Etot, evolution_T
del E_mean, E_fluct
del T_mean, T_fluct
del dummy
        
#%% EXERCISE 1.1 and 1.2 ----- Graphical stuff

# temp = np.loadtxt(folder_path + "/1_results_noPOLY_1000K.txt")
# dt, noPOLY_E, noPOLY_T = temp[:, 0], temp[:, 1:3], temp[:, 3:]

temp = np.loadtxt(folder_path + "/1_results_POLY_1000K.txt")
dt, POLY_E, POLY_T = temp[:, 0], temp[:, 1:3], temp[:, 3:]

temp = np.loadtxt(folder_path + "/1_results_POLY_2000K.txt")
dt, noPOLY_E, noPOLY_T = temp[:, 0], temp[:, 1:3], temp[:, 3:]

fig, ax = plt.subplots(1,1, figsize=(7,5.2))

ax.plot(dt, noPOLY_E[:, 1], "rx", label="2000 [K]")
ax.plot(dt, POLY_E[:,1], "bx", label="1000 [K]")

ax.set_xlabel("Time-step [fs]")
# ax.set_ylabel(r"$\sigma$(T) / $\langle T \rangle$")
ax.set_ylabel(r"$\sigma$(E) / $\langle E \rangle$")
# ax.set_title("Energy fluctuations vs Time-step size\n$T_{ini}$ = 1000 K")
ax.set_title("Energy fluctuations vs Time-step size\nPoly junction")
ax.grid()
ax.set_facecolor("whitesmoke")
ax.legend()
ax.set_yscale("log")

del temp
del dt
# del noPOLY_E, noPOLY_T
del POLY_E, POLY_T

#%% EXERCISE 1.3 ----- ANALYZE DATA

# sim_names = [
#     "50K_19fs_10ps",
#     "500K_5fs_10ps",
#     "1000K_3fs_10ps",
#     "1100K_2fs_10ps",
#     "1200K_2fs_10ps",
#     "1300K_2fs_10ps",
#     "1400K_2fs_10ps",
#     "1500K_2fs_10ps",
#     "1600K_2fs_10ps",
#     "1700K_2fs_10ps",
#     "1800K_2fs_10ps",
#     "1900K_2fs_10ps",
#     "2000K_2fs_10ps"]

# sim_T_ini = [50, 500, 1000, 1100, 1200, 1300, 1400, 1500,
#              1600, 1700, 1800, 1900, 2000]

sim_names = ["50K_2fs_30ps_PBC",   # <--- PBC checks
             "500K_2fs_30ps_PBC",
             "1000K_2fs_30ps_PBC",
             "2000K_2fs_30ps_PBC"]

sim_T_ini = [50, 500, 1000, 2000]   # <--- PBC checks

thermal_time = 10 * 1e-12

with open(folder_path + "/1_results_30ps_10ps_PBC.txt", "w") as f_results:

    for sim_name, T_ini in zip(sim_names, sim_T_ini):
        output_folder = folder_path + "/" + sim_name
        
        try:
            data = np.loadtxt(output_folder + "/Output.txt")
        except FileNotFoundError:
            del thermal_time
            del sim_names, sim_T_ini, sim_name, T_ini
            del output_folder, f_results
            del data
            del time, evolution_Etot, evolution_T
            del E_mean, E_sigma, E_fluct
            del T_mean, T_sigma, T_fluct
            raise
    
        time, evolution_Etot, evolution_T = data[:,1], data[:, 2], data[:, 5]

        E_mean, E_sigma, E_fluct = md.compute_RelFluct(evolution_Etot, time, verbose=False, thermal_time = thermal_time)
        T_mean, T_sigma, T_fluct = md.compute_RelFluct(evolution_T, time, verbose=False, thermal_time = thermal_time)
        
        print(sim_name)
        f_results.write(f"{str(T_ini).zfill(2)}\t{E_mean}\t{E_sigma}\t{T_mean}\t{T_sigma}\n")
        
del thermal_time
del sim_names, sim_T_ini, sim_name, T_ini
del output_folder, f_results
del data
del time, evolution_Etot, evolution_T
del E_mean, E_sigma, E_fluct
del T_mean, T_sigma, T_fluct

#%% EXERCISE 1.3 ----- Graphical stuff

fig, ax = plt.subplots(1,1, figsize=(7.3,6))

temp = np.loadtxt(folder_path + "/1_results_10ps_3ps.txt")
T_ini, T_eq, sigma_T = temp[:, 0], temp[:, 3], temp[:, 4]
ax.plot(T_ini, T_eq, "bx", label = r"10 ps - $t_{th}$ = 3 ps")
ax.errorbar(T_ini, T_eq, sigma_T, capsize = 5, fmt = ",b")

temp = np.loadtxt(folder_path + "/1_results_30ps_10ps.txt")
T_ini, T_eq, sigma_T = temp[:, 0], temp[:, 3], temp[:, 4]
ax.plot(T_ini, T_eq, "rx", label = r"30 ps - $t_{th}$ = 10 ps")
ax.errorbar(T_ini, T_eq, sigma_T, capsize = 5, fmt = ",r")

# temp = np.loadtxt(folder_path + "/1_results_30ps_10ps_PBC.txt")
# T_ini, T_eq, sigma_T = temp[:, 0], temp[:, 3], temp[:, 4]
# plt.plot(T_ini, T_eq, "ro", label = r"30 ps - $t_{th}$ = 10 ps - PBC")
# plt.errorbar(T_ini, T_eq, sigma_T, capsize = 5, fmt = ",r")


x_axis = np.linspace(0, 2000, 2001)
ax.plot(x_axis, 0.5*x_axis, "--k", label = "Expected", alpha = 0.5)

del temp
del T_ini
del T_eq, sigma_T
del x_axis

ax.set_xlabel(r"$T_{ini}$ [K]")
ax.set_ylabel(r"$T_{eq}$ [K]")
ax.set_title("Temperature behaviour\nEquilibrium vs Initial")
ax.grid()
ax.set_facecolor("whitesmoke")
ax.legend()
ax.set_yscale("linear")

#%% EXERCISE 1.4 and 1.5 ----- ANALYZE DATA

# ----- Convert the Conf files into Conf files with PBC for the adatom trajectory

# md.toPBC_ConfFile(old_folder = "point5_1300K_2fs_300ps")

# ----- Extract solely the adatom trajectory and save it in "Adatom_traj.txt"

# obtained_traj_file = md.extract_AdatomTraj(coords_folder = "point4_1700K_2fs_300ps")
obtained_traj_file = md.extract_AdatomTraj(coords_folder = "point4_1700K_PBC")
# obtained_traj_file = md.extract_AdatomTraj(coords_folder = "point5_1300K_2fs_300ps")
obtained_traj_file = md.extract_AdatomTraj(coords_folder = "point5_1300K_PBC")

#%% EXERCISE 1.4 and 1.5 ----- Graphical stuff

# ----- GRAFICI

traj_file = folder_path + "/point4_1700K_2fs_300ps/Adatom_traj.txt"
traj_file_PBC = folder_path + "/point4_1700K_PBC/Adatom_traj.txt"
# traj_file = folder_path + "/point5_1300K_2fs_300ps/Adatom_traj.txt"
# traj_file_PBC = folder_path + "/point5_1300K_PBC/Adatom_traj.txt"

traj = np.loadtxt(traj_file)
traj_PBC = np.loadtxt(traj_file_PBC)

# Color mapping

time = np.linspace(1, 1500, 1500)
norm = mpl.colors.Normalize(vmin=0, vmax=1500)
norm_time = norm(time)

cmap = mpl.colormaps['plasma']
colors = cmap(norm_time)

# Drawing

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6), width_ratios = [0.65, 0.35])
fig.suptitle("Adatom trajectory - FCC [1 0 0]", fontsize = 18, fontweight = 'bold')

ax1.scatter(traj[:, 0], traj[:, 1], c = colors, s = 0.9)

L = 20.3817
L = 16.6416
size = np.ones(1000)
scale = np.linspace(0, L, 1000)
ax1.scatter(size*scale, size*0, c = 'g', alpha = 0.1, s = 0.3)
ax1.scatter(size*L, size*scale, c = 'g', alpha = 0.05, s = 0.3)
ax1.scatter(size*scale, size*L, c = 'g', alpha = 0.05, s = 0.3)
ax1.scatter(size*0, size*scale, c = 'g', alpha = 0.1, s = 0.3)

ax1.set_title("Unwrapped trajectory")
ax1.set_xlabel("x-coord [A]")
ax1.set_ylabel("y-coord [A]")
ax1.set_aspect('equal')
ax1.set_facecolor('whitesmoke')
ax1.grid(True)

ax2.scatter(traj_PBC[:,0], traj_PBC[:,1], c = colors, s = 0.8)

ax2.set_title("Inside the sim cell")
ax2.set_xlabel("x-coord [A]")
ax2.set_ylabel("y-coord [A]")
ax2.set_ylim(-1, 17)
ax2.set_xlim(-1, 17)
ax2.set_aspect('equal')
ax2.set_facecolor('whitesmoke')
ax2.grid(True)

#%% EXERCISE 1.5 ----- Graphical stuff

dist = np.linspace(2, 5, 1000)
energy = md.LJ_6_12_energy(dist)

fig, ax = plt.subplots(1, 1)

ax.plot(dist, energy, 'b')

ax.vlines(4.5, -0.6, 2.1, ls = '--', color = 'k', alpha = 0.5)
ax.vlines(4.2, -0.6, 2.1, ls = '--', color = 'k', alpha = 0.5)

# ax.vlines(2.5477081738692124 / 2, -0.6, 2.1, ls = '--', color = 'k', alpha = 0.5)
# ax.vlines(4.1604 / 2, -0.6, 2.1, ls = '--', color = 'k', alpha = 0.5)

ax.set_title("Lennard-Jones 6-12 potential")
ax.set_xlabel("Distance [A]")
ax.set_ylabel("Energy [eV]")
ax.grid(True)
ax.set_facecolor("whitesmoke")
ax.set_ylim(-0.5, 2)

ax.annotate("Cut-off", xy = (4.52, 1.8))
ax.annotate("Poly", xy = (4.22, 1.8))
# ax.annotate(f"[1 0 0] basin\n{2.5477081738692124 / 2 :.2g}")
# ax.annotate("[1 1 1] basin\n{4.1604 / 2 : .2g}")

# Distance between two atoms in a fcc111 surface
# pos = np.array([[8.49238, 2.94185, 12.01],
#        [8.49238, 5.88369, 12.01]])

# i = int(0)
# j = int(1)

# dx = pos[i, 0] - pos[j, 0]
# dy = pos[i, 1] - pos[j, 1]
# dz = pos[i, 2] - pos[j, 2]
# dist2 = dx*dx + dy*dy + dz*dz
# dist = dist2**0.5

# print(dist)