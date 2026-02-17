# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 12:29:53 2026

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/MMC/Ag_surface_equilibrium/"

import numpy as np
import time
import os
from numba import jit

# ----- GLOBAL VARIABLES

LX, LY = 60, 60   # Surface size in atoms
N_ADATOMS = 25

J0, J1 = 4 * 0.345, 0.345   # Bonding energies in eV
K_BOLTZ = 1/11603   # In eV/K
TEMP = 100.   # Temperature in K

N_ITER = int(3 * 1e6)   # Number of iterations

SEED = 312428270054674116253981409899306192348   # Obtained with secrets.randbits(128)
TEMP_BOOL = False
OUTPUT_FOLDER = ""   # Folder storing the entirety of the sim's results


#%% -------------------- FUNCTION DEFINITIONS -------------------- %%#

# ----- I/O function

# The 'Input_parameters.txt' file must be composed as:
    # 'TEMP_BOOL' a bool deciding whether to use Metropolis selection (T !=0) or not (T == 0).
    #                  Note that it also selects the proper 'OUTPUT_FOLDER'.
    # 'TEMP' value in Kelvin
    # 'N_ADATOMS' the number of moving adatoms on the surface
    # 'SEED' value used for random extractions, must be an INTEGER
    # 'output_name' name of the folder where the output is sent
    
    # 'TEMPORARY' used differently in different simulations
    
# Output values:
    # 'f_output' a file object, which is open and MUST BE CLOSED manually at the end of the simulation

def read_input (input_file: str = "Input_parameters.txt"):
    
    # ----- Read the input parameters from the 'input_file' ------
    
    filename = full_path + input_file
    
    with open(filename, 'r') as input_params:
    
        global OUTPUT_FOLDER
        global N_ADATOMS
        global TEMP_BOOL
        global TEMP
        global SEED
    
        TEMP_BOOL = input_params.readline()
        TEMP_BOOL = TEMP_BOOL.splitlines()[0]   # Remove the '\n' at the end
        
        TEMP = float(input_params.readline())
        N_ADATOMS = int(input_params.readline())
        SEED = int(input_params.readline())
        
        if TEMP_BOOL == "True":
            TEMP_BOOL = True
            OUTPUT_FOLDER = "2nd_part"
        else:
            TEMP_BOOL = False
            OUTPUT_FOLDER = "1st_part"
            
        output_name = input_params.readline()
        output_name = output_name.splitlines()[0]
            
        # time_stamp = time.strftime("%m_%b%d_%H_%M", time.localtime())   # Time format: month_day_hour_minute
        # OUTPUT_FOLDER = full_path + OUTPUT_FOLDER + "/" + time_stamp + "_" + output_name
        OUTPUT_FOLDER = full_path + OUTPUT_FOLDER + "/" + output_name
        
        try:
            os.mkdir(OUTPUT_FOLDER)
        except FileExistsError:
            if (input("The output folder already exists. Rewrite it [y/N]?\t") == "N"):
                print("\n")
                raise 
                
        f_output = open(OUTPUT_FOLDER + "/Output.txt", 'w')
        
        # global LX
        # global LY
        # global TEMPORARY
        
        # TEMPORARY = float(input_params.readline())
        # LX = TEMPORARY
        # LY = TEMPORARY
        
    return f_output

# ----- Save a snap of the actual atomic configuration. Compatible with OVITO.

def save_snap (filename: str, i: int, height: np.array) -> None:
    
    global OUTPUT_FOLDER
    
    with open(OUTPUT_FOLDER + "/" + filename, 'w') as f_snap:
        f_snap.write(f"{np.sum(height+1)} \n{i} \n")
        
        for x in range(0, LX):
            for y in range(0, LY):
                for z in range(0, height[x, y] + 1):
                    f_snap.write(f"{x} \t {y} \t {z} \n")
                    
    return None
    
# ----- Useful routines:
    # 1. Computes the energy of the entire surface configuration
    # 2. Converts a generic atomic position into one allowable by PBCs.
    # 3. Counts the number of nbrs around the atom (x, y)
    # 4. Updates the 'nnbrs' matrix after a deposition event in (xdepo, ydepo)
    # 5. Updates the 'nnbrs' matrix after a diffusion event (x0, y0) -> (xdiff, ydiff)

@jit   # FUNCTION 1
def compute_energy (height: np.array, nnbrs: np.array) -> float:

    global LX
    global LY
    global J0
    global J1
    
    energy = 0
    
    for x in range(0, LX):
        for y in range(0, LY):
            for z in range(1, height[x, y] + 1):
                energy += - J0 - 0.5 * nnbrs[x, y] * J1
    return energy

@jit   # FUNCTION 2
def into_PBC (coord: int, L: int) -> int:
    
    coord = int(coord - L * np.floor(coord/L))
    
    return coord

@jit   # FUNCTION 3
def count_nbrs(x: int, y: int, height: np.array) -> int:
    
    n = 0
    z = height[x, y]
    
    global LX
    global LY
    
    # Coordinates of the atoms around (x, y) converted with the PBCs
    
    x_east, x_west = into_PBC(x + 1, LX), into_PBC(x - 1, LX)
    y_south, y_north = into_PBC(y - 1, LY), into_PBC(y + 1, LY)
        
    if (z <= height[x_east, y]):   # East
        n += 1
    if (z <= height[x, y_south]):   # South
        n += 1
    if (z <= height[x_west, y]):   # West
        n += 1
    if (z <= height[x, y_north]):   # North
        n += 1
        
    return n

@jit   # FUNCTION 4
def update_depo_nnbrs (xdepo: int, ydepo: int, height: np.array, nnbrs: np.array) -> np.array:
    
    global LX
    global LY
    
    # Coordinates of the atoms around (xdepo, ydepo) converted with the PBCs 
    
    x_east, x_west = into_PBC(xdepo + 1, LX), into_PBC(xdepo - 1, LX)
    y_south, y_north = into_PBC(ydepo - 1, LY), into_PBC(ydepo + 1, LY)
    
    nnbrs[xdepo, ydepo] = count_nbrs(xdepo, ydepo, height)   # Deposition location
    nnbrs[x_east, ydepo] = count_nbrs(x_east, ydepo, height)   # East
    nnbrs[x_west, ydepo] = count_nbrs(x_west, ydepo, height)   # West
    nnbrs[xdepo, y_south] = count_nbrs(xdepo, y_south, height)   # South
    nnbrs[xdepo, y_north] = count_nbrs(xdepo, y_north, height)   # North
    
    return nnbrs

@jit   # FUNCTION 5
def update_nnbrs (x0: int, y0: int, xdiff: int, ydiff: int, 
                       height: np.array, nnbrs: np.array) -> np.array:
    
    # Updates each position using the related deposition function.
    # It is somehow redundant because (x0, y0) and (xdiff, ydiff) may be updated twice,
    # but I am lazy in working out a proper and slightly faster solution...
    
    nnbrs = update_depo_nnbrs(x0, y0, height, nnbrs)
    nnbrs = update_depo_nnbrs(xdiff, ydiff, height, nnbrs)
    
    return nnbrs

#%% -------------------- MAIN CODE -------------------- %%#

# ----- FUNDAMENTAL VARIABLES: initializations

f_output = read_input()

RNG = np.random.default_rng(seed=SEED)
height = np.zeros((LX, LY), dtype = np.int_)
nnbrs = np.zeros((LX, LY), dtype = np.int_)

# ----- DEPOSITION EVENTS

for i in range(0, N_ADATOMS):
    xdepo = RNG.integers(0, LX, endpoint=False)
    ydepo = RNG.integers(0, LY, endpoint=False)
    height[xdepo, ydepo] += 1
    nnbrs = update_depo_nnbrs(xdepo, ydepo, height, nnbrs)
    
old_Energy = compute_energy(height, nnbrs)
f_output.write(f"{0:d} \t\t {old_Energy:.21E} \t {1:d}\n")   # Save the initial configuration as well

# NOTE: 'height' cannot be < 0, because the algorithm chooses to move only from columns
#       where 'height' is strictly > 0.

# ----- OUTPUT RELATED VARIABLES

snap_di = int(50 * 1e3)
snap_cumulant = 0
n_snap = 0

is_accepted = 1   # Keep track of the acceptance: 1 == trial accepted, 0 == trial rejected

# ----- MC LOOP: T = 0 K

if (TEMP_BOOL == False):
    print("\nT == 0 detected\n\n")
    
    for i in range(1, N_ITER + 1):
                
        # ----- SELECTION OF THE EVENT
        
        coord_list = np.argwhere(height > 0)   # A list of (x, y) coords
        chosen_coord = RNG.integers(0, len(coord_list), endpoint=False)
        xold, yold = coord_list[chosen_coord]
        
        xnew = RNG.integers(0, LX, endpoint=False)
        ynew = RNG.integers(0, LY, endpoint=False)
        
        # ----- TRIAL MOVE
        
        height[xold, yold] -= 1
        height[xnew, ynew] += 1
        
        nnbrs = update_nnbrs(xold, yold, xnew, ynew, height, nnbrs)
        
        new_Energy = compute_energy(height, nnbrs)
        
        # ----- ACCEPTANCE OF THE EVENT
        dE = new_Energy - old_Energy
        
        if (dE < 1e-8):
            old_Energy = new_Energy
            
            is_accepted = 1
            
        else:   # Restore the previous configuration
            height[xold, yold] += 1
            height[xnew, ynew] -= 1
            
            nnbrs = update_nnbrs(xold, yold, xnew, ynew, height, nnbrs)
            
            is_accepted = 0
        
        # ----- END OF STEP ROUTINES: compute/save outputs
        
        # Regular output (at each step)
        
        f_output.write(f"{i:d} \t\t {old_Energy:.21E} \t {is_accepted:d}\n")
        
        # Save a snap of the atomic configuration every 'snap_dt'
        
        if (i > snap_cumulant):
            save_snap("SNAP" + str(n_snap).zfill(4) + ".xyz", i, height)
            
            n_snap += 1
            snap_cumulant += snap_di
            
            print(f"{i:E} / {N_ITER:E}")
        
# ----- MC LOOP: T != 0 K

elif (TEMP_BOOL == True):
    print("\nT > 0 detected\n\n")
    
    for i in range(1, N_ITER + 1):
                
        # ----- SELECTION OF THE EVENT
        
        coord_list = np.argwhere(height > 0)   # A list of (x, y) coords
        chosen_coord = RNG.integers(0, len(coord_list), endpoint=False)
        xold, yold = coord_list[chosen_coord]
        
        xnew = RNG.integers(0, LX, endpoint=False)
        ynew = RNG.integers(0, LY, endpoint=False)
        
        # ----- TRIAL MOVE
        
        height[xold, yold] -= 1
        height[xnew, ynew] += 1
        
        nnbrs = update_nnbrs(xold, yold, xnew, ynew, height, nnbrs)
        
        new_Energy = compute_energy(height, nnbrs)
        
        # ----- ACCEPTANCE OF THE EVENT
        dE = new_Energy - old_Energy
        
        if (dE < 1e-8):
            old_Energy = new_Energy
            
            is_accepted = 1
            
        elif (RNG.random() < np.exp(- dE /(K_BOLTZ*TEMP) ) ):
            old_Energy = new_Energy
            
            is_accepted = 1
            
        else:   # Restore the previous configuration
            height[xold, yold] += 1
            height[xnew, ynew] -= 1
            
            nnbrs = update_nnbrs(xold, yold, xnew, ynew, height, nnbrs)
            
            is_accepted = 0
        
        # ----- END OF STEP ROUTINES: compute/save outputs
        
        # Regular output (at each step)
        
        f_output.write(f"{i:d} \t\t {old_Energy:.21E} \t {is_accepted:d}\n")
        
        # Save a snap of the atomic configuration every 'snap_dt'
        
        if (i > snap_cumulant):
            # save_snap("SNAP" + str(n_snap).zfill(4) + ".xyz", i, height)
            
            n_snap += 1
            snap_cumulant += snap_di
            
            print(f"{i:E} / {N_ITER:E}")
    
# ----- END OF SIMULATION ROUTINES: things to do before ending code execution

f_output.close()

with open(OUTPUT_FOLDER + "/YOURinput_parameters.txt", "w") as input_params:
    input_params.write(str(TEMP_BOOL) + "\n")
    input_params.write(str(TEMP) + "\n")
    input_params.write(str(N_ADATOMS) + "\n")
    input_params.write(str(SEED) + "\n")
    
#%%
f_output.close()
