# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 14:43:15 2025

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/"
# folder names: "Working_folder/"   "Useful_files_for_simulations/"

import sys
sys.path.append(full_path)   # <- necessary when using 'launch_several_sims.py'
import my_package.my_module as md

import numpy as np
import matplotlib.pyplot as plt
import os

#%% MAIN CODE

# Things that you can change IN THE CODE (not in 'Input_parameters'):
#   - # of nbrs in 'whichnbr'
#   - nPrint

# --------------- INITIALIZATIONS --------------

pos, Temp, myseed, output_folder = md.read_input("Input_parameters.txt")
natoms = int(np.shape(pos)[0])
whichnbr = np.zeros((natoms, 30), dtype=int)   # nbrs matrix of each atom

# NOTE: an fcc lattice has at most 12 nn. 
#       If 'rc' includes dnn2, then 12 MUST be suitably changed to 30 (17 would've been enough) in 'whichnbr'.

# OUTPUT and PRINTING infos
f_output = open(output_folder + "/Output.txt", 'w')   # <--- Needed for post-processing
# f_status = open(output_folder + "/Status.txt", 'a')
nPrint = 100   # <--- Print simulation infos every 'nPrint' timesteps
nOut = 0   # Dummy variable for saving lattice's configs

dt = md.TEMPORARY * 1e-15
n_timesteps = int( np.ceil(300 * 1e-12 / dt) )

# --------------- COMPUTATIONS ---------------

# nnbrs: lists the number of nbrs of each atom (array of dimension 'natoms')
# force: lists the forces acting on each atom (matrix of dimension 'natoms' x 3)
# Epot: total potential energy of the system
# vel: lists the velocities of each atom (matrix of dimension 'natoms' x 3)
# Ekin: total kinetic energy of the system

nnbrs, whichnbr = md.compute_fastVerlet(whichnbr, pos, natoms)
force = md.compute_force(nnbrs, whichnbr, pos, natoms)
vel = md.compute_initialVelocity(Temp, natoms, myseed)

for iStep in range(0, n_timesteps):
    
    time = iStep * dt
    
    # ----- UPDATE POSITIONS and VELOCITIES -----
    
    pos = pos + vel*dt + 0.5 * (force/md.MASS) * dt**2
    
    nnbrs, whichnbr = md.compute_fastVerlet(whichnbr, pos, natoms)
    force_next = md.compute_force(nnbrs, whichnbr, pos, natoms)
    
    vel = vel + 0.5 * ( (force + force_next) / md.MASS ) * dt
    
    force = force_next
    
    # ----- PHYSICAL VALUES -----
    
    Ekin = md.compute_fastEkin(vel)
    Epot = md.compute_Epot(nnbrs, whichnbr, pos, natoms)
    Etot = Ekin + Epot
    
    inst_Temp = Ekin / (natoms * md.K_BOLTZ * 1.5)
    
    f_output.write(str(iStep).zfill(4) + f" \t {time:.4e} \t {Etot} \t {Epot} \t {Ekin} \t {inst_Temp}\n")
    
    # Save the lattice configuration and output the simulation status
    if ( (iStep % nPrint) == (nPrint - 1)):
        nOut += 1
        coords_file = output_folder + "/Conf" + str(nOut).zfill(3) + ".xyz"
        np.savetxt(coords_file, pos, header = f"{natoms}\nTime={time}", comments="")
        
        # Use THIS when launching from Spyder Console or
        # when launching with a secondary script from Anaconda Powershell
        print(f"iStep = {iStep:n} \t\t Time = {time:.2E} \t Energy = {Etot}\n")
        
        # Use THESE when launching with a secondary script from Spyder Console
        # f_status.write(f"iStep = {iStep:n} \t\t Time = {time:.2E} \t Energy = {Etot}\n")
        # f_status.flush()
    
np.savetxt(output_folder + "/YOURInput_parameters.txt", [len(pos), md.RC, md.RP, Temp, 
                                                         md.L_PBC[0], md.L_PBC[1], md.L_PBC[2], myseed])
    
f_output.close()
# f_status.close()

# os.remove(output_folder + "/Status.txt")

#%% TEST

