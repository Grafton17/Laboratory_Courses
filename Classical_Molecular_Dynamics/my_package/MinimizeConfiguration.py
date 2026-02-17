# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 22:50:26 2025

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/"   # Where the module sits

import numpy as np
import os
import matplotlib.pyplot as plt
import my_package.my_module as md

C_steep = 0.005
limForce = 1 * 1e-4   # Typically at 1 * 1e-3, set to 1e-4 for 'fcc111a336+1.txt'
n_timesteps = 5000   # Typically at 3000, set to 5000 for 'fcc111a336+1.txt'

#%%

# The structure of "Input_parameters.txt" is illustrated in "my_module.py":
    # give the proper 'LATTICE_FILE'
    # set '0 0 0' in 'use_PBCs' as we have to minimize the cluster configuration
    # set 'use_junction' to True
    # set a proper 'output_folder' name
    # all other parameters are not important

# In this context, 'Temp' and 'myseed' from md.read_input() are not used

all_maxForce = []
all_indexes = []
all_Epot = []

pos, Temp, myseed, output_folder = md.read_input("Input_parameters.txt")
natoms = int(np.shape(pos)[0])
whichnbr = np.zeros((natoms, 30), dtype=int)   # nbrs matrix of each atom

for iStep in range(0, n_timesteps):
    nnbrs, whichnbr = md.compute_fastVerlet(whichnbr, pos, natoms)
    force = md.compute_force(nnbrs, whichnbr, pos, natoms)
    
    maxForce = np.max( np.sqrt( force[:, 0]**2 + force[:, 1]**2 + force[:, 2]**2 ) )
    ind_maxForce = np.argmax( np.sqrt( force[:, 0]**2 + force[:, 1]**2 + force[:, 2]**2 ) )
    Epot = md.compute_Epot(nnbrs, whichnbr, pos, natoms)
    
    all_maxForce.append(maxForce)
    all_indexes.append(ind_maxForce)
    all_Epot.append(Epot)
    
    print(f"{iStep:n}\t\t{ind_maxForce:n}\t\t{maxForce:.2E}\t{Epot}", flush = True)
    
    if maxForce < limForce:
        print("\nEnd")
        break
    else:
        pos = pos + C_steep * force
    
np.savetxt(output_folder + f"/MINIMIZED_" + md.LATTICE_FILE, pos,
           header = f"{natoms}\n", comments = "")
np.savetxt(output_folder + f"/MINIMIZED_infos_" + md.LATTICE_FILE, 
           np.stack((all_maxForce, all_indexes, all_Epot), axis=-1),
           header = f"Csteep = {C_steep:.2E}\tlimForce = {limForce:.2E}\nmaxForce - maxForce atom - Epot", 
           comments = "", fmt = '%.18e, %5u, %.18e')
    
del pos
del Temp
del myseed
del natoms
del whichnbr
del nnbrs
del force
del maxForce
del ind_maxForce
del Epot

#%% Plot and save energy and maxForce vs iteration step in the same output folder

fig, (axE, axF) = plt.subplots(1, 2, figsize=(10,5), layout='constrained')
# fig.tight_layout(h_pad = 5)

x_range = range(0, iStep+1)

axE.plot(x_range, all_Epot, '-')
axE.set_title("Total energy")
axE.set_ylabel("Energy [eV]")
axE.grid()
axE.ticklabel_format(axis='y', useOffset = all_Epot[0], style = 'plain', useMathText = True)

axF.plot(x_range, all_maxForce, '-')
axF.set_title("Max force")
axF.set_ylabel("Force [eV/A]")
axF.grid()

fig.supxlabel("# of iteration")

plt.savefig(output_folder + f"/MINIMIZED_" + md.LATTICE_FILE + ".pdf")
# plt.show()
plt.close()

#%%
a = np.stack((all_maxForce, all_indexes, all_Epot), axis=-1)