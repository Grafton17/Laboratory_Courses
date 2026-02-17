# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 12:30:13 2026

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/MMC/Ag_surface_equilibrium/"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.constants as cst

import os
import sys

mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.facecolor'] = 'whitesmoke'
mpl.rcParams['figure.titlesize'] = 20
mpl.rcParams['figure.titleweight'] = 'bold'
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['figure.figsize'] = (8, 6)
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.bbox'] = 'tight'

#%% FUNCTIONS

def compute_RelFluct (data: np.ndarray, time: np.ndarray, thermal_time: int = 5 * 1e5, 
                      verbose = True) -> (float, float):
    """Compute the mean, absolute and relative fluctuation of 'data' after mixing has occured."""
    
    thermal_index = 0

    for i in range(0, len(data)):
        if (time[i] > thermal_time ):
            thermal_index = i
            break
        
    data_mean = np.mean(data[thermal_index:])

    data2_mean = np.mean(data[thermal_index:]**2)
    fluct = np.sqrt(data2_mean - data_mean**2)
    # E_fluct = np.std(evolution_Etot, mean=E_mean)   # <--- Available in newer Numpy versions

    rel_fluct = fluct / np.abs(data_mean)
    
    if verbose:
        print(f"<S>: {data_mean} \t Std: {fluct} \t rel_fluct: {rel_fluct:.3E}")
    
    return data_mean, fluct, rel_fluct

#%% POINT 1 ----- Data analysis

sim_names = ["2_0_60_3", "3_0_60_3", "4_0_60_3", "5_0_60_3", "6_0_60_3", "7_0_60_3", 
              "8_0_60_3", "9_0_60_3", "10_0_60_3", "11_0_60_3", "12_0_60_3", "13_0_60_3", 
              "14_0_60_3", "15_0_60_3", "16_0_60_3", "17_0_60_3", "18_0_60_3", "19_0_60_3", 
              "20_0_60_3"]

sim_N = [ 2,  3,  4,  5,  6,  7,  
          8, 9 , 10, 11, 12, 13, 
          14, 15, 16, 17, 18, 19, 
          20]

means = np.array([])
stds = np.array([])

for folder_name in sim_names:

    folder_path = full_path + "1st_part/" + folder_name + "/Output.txt"

    data = np.loadtxt(folder_path)
    step, energy = data[:,0], data[:,1]
    
    mean, std, fuffa = compute_RelFluct(energy, step, thermal_time = 1 * 1e5, verbose = True)
    means = np.append(means, mean)
    stds = np.append(stds, std)
    
    print(folder_name)

del data
del step
del energy

save_file = full_path + "1st_part/Energy_Chem_vs_N_seed3.txt"

np.savetxt(save_file, np.stack((means, means/sim_N, sim_N), axis = -1))

#%% POINT 1 ----- Graphical stuff

# ------ SINGLE SIM CHECK

# folder_name = "4_0_60_1"
# folder_path = full_path + "1st_part/" + folder_name + "/Output.txt"

# cut = int(1e6)

# data = np.loadtxt(folder_path)
# step, energy = data[:cut,0], data[:cut,1]

# fig, ax = plt.subplots(1, 1, figsize=(8,7))

# ax.plot(step, energy)

# ax.grid(True)
# ax.set_facecolor("whitesmoke")
# ax.set_ylabel("Energy [eV]")
# ax.set_xlabel("MC step")
# ax.set_title("Energy evolution: seed #1\nN = 4 - T = 0 [K]")
# ax.set_xlim(-1500, 7000)

# ----- ENERGY and CHEMICAL POTENTIAL

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6.3), sharex = True)

ax1.set_title(r"Comparison: $E_{min} (N)$ vs $N$")
# ax1.set_xlabel("Number of adatoms")
ax1.set_ylabel(r"$E_{min} (N)$ [eV]")

ax2.set_title(r"Comparison: $\mu (N)$ vs $N$")
ax2.set_xlabel("Number of adatoms")
ax2.set_ylabel(r"$\mu (N)$ [eV]")
ax2.ticklabel_format(axis = 'x', style = 'plain')

for i, c in zip([1, 3], ['royalblue', 'darkred']):

    save_file = full_path + f"1st_part/Energy_Chem_vs_N_seed{i:d}.txt"
    data = np.loadtxt(save_file)
    energy, chem, Ns = data[:, 0], data[:, 1], data[:, 2]
    
    ax1.plot(Ns, energy, 'o--', c = c, label = f"Seed #{i:d}")
    
    ax2.plot(Ns, chem, 'o--', c = c, label = f"Seed #{i:d}")

ax1.legend()
ax2.legend()
fig.show()

#%% POINT 2 ----- Data analysis

# sim_names = ["25_1_60_1", "25_10_60_1", "25_100_60_1", "25_500_60_1", "25_1000_60_1", "25_1500_60_1", 
#              "25_2000_60_1", "25_2500_60_1", "25_2700_60_1"]

# sim_T = [1, 10, 100, 500, 1000, 1500,
#          2000, 2500, 2700]

# sim_names = ["25_1_60_1", "25_10_60_1", "25_100_60_1", "25_200_60_1", "25_300_60_1", "25_500_60_1", 
#               "25_600_60_1", "25_750_60_1", "25_900_60_1", "25_1000_60_1", "25_1100_60_1", "25_1200_60_1",
#               "25_1300_60_1", "25_1400_60_1", "25_1500_60_1", "25_1600_60_1", "25_1750_60_1", "25_1900_60_1",
#               "25_2000_60_1", "25_2250_60_1", "25_2500_60_1", "25_2700_60_1"]

# sim_names = ["25_1_60_3", "25_10_60_3", "25_100_60_3", "25_200_60_3", "25_300_60_3", "25_500_60_3", 
#               "25_600_60_3", "25_750_60_3", "25_900_60_3", "25_1000_60_3", "25_1100_60_3", "25_1200_60_3", "25_1250_60_3",
#               "25_1300_60_3", "25_1400_60_3", "25_1500_60_3", "25_1600_60_3", "25_1750_60_3", "25_1900_60_3",
#               "25_2000_60_3", "25_2250_60_3", "25_2500_60_3", "25_2700_60_3"]

sim_names = ["25_1_40_3", "25_10_40_3", "25_100_40_3", "25_200_40_3", "25_300_40_3", "25_500_40_3", 
              "25_600_40_3", "25_750_40_3", "25_900_40_3", "25_1000_40_3", "25_1100_40_3", "25_1200_40_3", "25_1250_40_3",
              "25_1300_40_3", "25_1400_40_3", "25_1500_40_3", "25_1600_40_3", "25_1750_40_3", "25_1900_40_3",
              "25_2000_40_3", "25_2250_40_3", "25_2500_40_3", "25_2700_40_3"]

sim_T = [1, 10, 100, 200, 300, 500, 
            600, 750, 900, 1000, 1100, 1200, 1250,
            1300, 1400, 1500, 1600, 1750, 1900, 
            2000, 2250, 2500, 2700]

# sim_names = ["25_1_20_3", "25_10_20_3", "25_100_20_3", "25_200_20_3", "25_300_20_3", "25_500_20_3", 
#               "25_600_20_3", "25_750_20_3", "25_900_20_3", "25_1000_20_3", "25_1100_20_3", "25_1200_20_3",
#               "25_1300_20_3", "25_1400_20_3", "25_1500_20_3", "25_1600_20_3", "25_1750_20_3", "25_1900_20_3",
#                 "25_2000_20_3", "25_2250_20_3", "25_2500_20_3", "25_2700_20_3"]


# sim_T = [1, 10, 100, 200, 300, 500, 
#             600, 750, 900, 1000, 1100, 1200, 
#             1300, 1400, 1500, 1600, 1750, 1900, 
#             2000, 2250, 2500, 2700]

# sim_names = ["25_1_60_5", "25_10_60_5", "25_100_60_5", "25_200_60_5", "25_300_60_5", "25_500_60_5", 
#               "25_600_60_5", "25_750_60_5", "25_900_60_5", "25_1000_60_5", "25_1100_60_5", "25_1200_60_5",
#               "25_1300_60_5", "25_1400_60_5", "25_1500_60_5", "25_1600_60_5", "25_1750_60_5", "25_1900_60_5",
#               "25_2000_60_5", "25_2250_60_5", "25_2500_60_5", "25_2700_60_5"]

# sim_T = [1, 10, 100, 200, 300, 500, 
#             600, 750, 900, 1000, 1100, 1200, 
#             1300, 1400, 1500, 1600, 1750, 1900, 
#             2000, 2250, 2500, 2700]

means = np.array([])
stds = np.array([])

for folder_name in sim_names:

    folder_path = full_path + "2nd_part/" + folder_name + "/Output.txt"

    data = np.loadtxt(folder_path)
    step, energy = data[:,0], data[:,1]
    
    mean, std, fuffa = compute_RelFluct(energy, step, thermal_time = 10 * 1e5, verbose = True)
    means = np.append(means, mean)
    stds = np.append(stds, std)
    
    print(folder_name)

del data
del step
del energy

save_file = full_path + "2nd_part/Energy_Std_vs_T_seed3_L40.txt"

np.savetxt(save_file, np.stack((means, stds, sim_T), axis = -1))

#%% POINT 2 ----- Graphical stuff

# ------ SINGLE SIM CHECK

# folder_name = "25_1500_60_1"
# folder_path = full_path + "2nd_part/" + folder_name + "/Output.txt"
# cut = int(3e6)

# data = np.loadtxt(folder_path)
# step, energy = data[:cut,0], data[:cut,1]

# fig, ax = plt.subplots(1, 1, figsize=(8,7))

# ax.plot(step, energy)

# ax.grid(True)
# ax.set_facecolor("whitesmoke")
# ax.set_ylabel("Energy [eV]")
# ax.set_xlabel("MC step")
# ax.set_title("Energy evolution: seed #1\nN = 25 - T = ... [K]")
# ax.set_xlim(-1500, 7000)






# ----- ENERGY and FLUCTUATIONS with/without seed comparison

# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex = True)

# ax1.set_title(r"Mean energy vs Temperature")
# # ax1.set_xlabel("Number of adatoms")
# ax1.set_ylabel(r"$\langle E \rangle$ [eV]")

# ax2.set_title(r"Fluctuations vs Temperature")
# ax2.set_xlabel("Temperature [K]")
# ax2.set_ylabel(r"$\sigma_{E}$ [eV]")
# ax2.ticklabel_format(axis = 'x', style = 'plain')

# fig.suptitle("60 x 60 lattice")

# for i, c in zip([1, 3, 5], ['royalblue', 'darkred', 'darkgreen']):
    
#     save_file = full_path + f"2nd_part/Energy_Std_vs_T_seed{i:d}_L60.txt"
#     data = np.loadtxt(save_file)
#     energy, std, Ts = data[:, 0], data[:, 1], data[:, 2]

#     ax1.plot(Ts, energy, 'o--', label = f"Seed #{i:d}", c = c)

#     ax2.plot(Ts, std, 'o--', label = f"Seed #{i:d}", c = c)

# ax1.legend()
# ax2.legend()
# fig.show()
    
    
    
# ----- SEVERAL SEEDS COMPARISON - evolution

# fig, ax = plt.subplots(1, 1, figsize=(8,7))

# for i in [1, 3, 5]:
#     folder_name = "25_100_60_{i:d}"
#     folder_path = full_path + "2nd_part/" + folder_name + "/Output.txt"
#     cut = int(1e6)

#     data = np.loadtxt(folder_path)
#     step, energy = data[:cut,0], data[:cut,1]

#     ax.plot(step, energy, label = f"Seed #{i:d}")

# ax.grid(True)
# ax.set_facecolor("whitesmoke")
# ax.set_ylabel("Energy [eV]")
# ax.set_xlabel("MC step")
# ax.set_title("Energy evolution: seed #3\nN = 25 - T = 500 [K]")
# ax.legend()
# ax.set_xlim(-1500, 7000)





# ----- ENERGY and FLUCTUATIONS with different lattices

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex = True)

ax1.set_title(r"Mean energy vs Temperature")
# ax1.set_xlabel("Number of adatoms")
ax1.set_ylabel(r"$\langle E \rangle$ [eV]")

ax2.set_title(r"Fluctuations vs Temperature")
ax2.set_xlabel("Temperature [K]")
ax2.set_ylabel(r"$\sigma_{E}$ [eV]")
ax2.ticklabel_format(axis = 'x', style = 'plain')

fig.suptitle("Seed #3")

for i, c in zip([60, 40, 20], ['royalblue', 'darkred', 'darkgreen']):
    
    save_file = full_path + f"2nd_part/Energy_Std_vs_T_seed3_L{i:d}.txt"
    data = np.loadtxt(save_file)
    energy, std, Ts = data[:, 0], data[:, 1], data[:, 2]

    ax1.plot(Ts, energy, 'o--', label = f"L = {i:d}", c = c)

    ax2.plot(Ts, std, 'o--', label = f"L = {i:d}", c = c)

ax1.legend()
ax2.legend()
fig.show()



#%% 

