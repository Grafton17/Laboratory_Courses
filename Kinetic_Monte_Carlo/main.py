# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 13:26:16 2026

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/KMC/Ag epitaxy/"

import numpy as np
import time
import os
from numba import jit

# ----- GLOBAL VARIABLES

LX, LY = 60, 60   # Surface size in atoms

J0, J1 = 4 * 0.345, 0.345   # Bonding energies in eV
K_BOLTZ = 1/11603   # In eV/K
NU = 1 * 1e13   # Attempt frequency in Hz
TEMP = 100.   # Temperature in K

PHI = 0.2   # Deposition rate in ML/s (ie '# of atoms per second in each column')s
THETA = 5.   # Nominal coverage (to reach) in ML/s
K_DEPO = PHI * LX * LY   # Deposition (transition) rate

END_TIME = THETA / PHI   # Time to simulate in s

SEED = 312428270054674116253981409899306192348   # Obtained with secrets.randbits(128)
DIFFUSION_BOOL = False
OUTPUT_FOLDER = ""   # Folder storing the entirety of the sim's results


#%% -------------------- FUNCTION DEFINITIONS -------------------- %%#

# ----- I/O function

# The 'Input_parameters.txt' file must be composed as:
    # 'DIFFUSION_BOOL' a bool deciding whether to activate diffusion or not.
    #                  Note that it also selects the proper 'OUTPUT_FOLDER'.
    # 'PHI' the deposition rate, value in ML/s
    # 'TEMP' value in Kelvin
    # 'SEED' value used for random extractions, must be an INTEGER
    
    # 'TEMPORARY' used differently in different simulations
    
# Output values:
    # 'f_output' a file object, which is open and MUST BE CLOSED manually at the end of the simulation

def read_input (input_file: str = "Input_parameters.txt"):
    
    # ----- Read the input parameters from the 'input_file' ------
    
    filename = full_path + input_file
    
    with open(filename, 'r') as input_params:
    
        global OUTPUT_FOLDER
        global DIFFUSION_BOOL
        global PHI
        global K_DEPO
        global TEMP
        global SEED
    
        DIFFUSION_BOOL = input_params.readline()
        DIFFUSION_BOOL = DIFFUSION_BOOL.splitlines()[0]   # Remove the '\n' at the end
        
        PHI = float(input_params.readline())
        K_DEPO = PHI * LX * LY
        TEMP = float(input_params.readline())
        SEED = int(input_params.readline())
        
        if DIFFUSION_BOOL == "True":
            DIFFUSION_BOOL = True
            OUTPUT_FOLDER = "2nd_part"
        else:
            DIFFUSION_BOOL = False
            OUTPUT_FOLDER = "1st_part"
            
        time_stamp = time.strftime("%m_%b%d_%H_%M", time.localtime())   # Time format: month_day_hour_minute
        OUTPUT_FOLDER = full_path + OUTPUT_FOLDER + "/" + time_stamp + f"_phi_{PHI}"
        
        try:
            os.mkdir(OUTPUT_FOLDER)
        except FileExistsError:
            if (input("The output folder already exists. Rewrite it [y/N]?\t") == "N"):
                print("\n")
                raise 
                
        f_output = open(OUTPUT_FOLDER + "/Output.txt", 'w')
        
        # TEMPORARY = float(input_params.readline())
        
    return f_output

# ----- Save a snap of the actual atomic configuration. Compatible with OVITO.

def save_snap (filename: str, snap_actual_coverage: float, height: np.array) -> None:
    
    global OUTPUT_FOLDER
    
    with open(OUTPUT_FOLDER + "/" + filename, 'w') as f_snap:
        f_snap.write(f"{np.sum(height+1)} \n{snap_actual_coverage} \n")
        
        for x in range(0, LX):
            for y in range(0, LY):
                for z in range(0, height[x, y] + 1):
                    f_snap.write(f"{x} \t {y} \t {z} \n")
                    
    # f_snap = open(OUTPUT_FOLDER + "/" + filename, 'w')
    # f_snap.write(f"{np.sum(height+1)} \n{snap_actual_coverage} \n")
    
    # for x in range(0, LX):
    #     for y in range(0, LY):
    #         for z in range(0, height[x, y] + 1):
    #             f_snap.write(f"{x} \t {y} \t {z} \n")
                
    # f_snap.close()
    
# ----- Useful routines:
    # 1. Given 'rho' and the 'K_DIFF' matrix it finds the chosen adatom.
    # 2. Converts a generic atomic position into one allowable by PBCs.
    # 3. Counts the number of nbrs around the atom (x, y)
    # 4. Updates the 'nnbrs' matrix after a deposition event in (xdepo, ydepo)
    # 5. Updates the 'nnbrs' matrix after a diffusion event (x0, y0) -> (xdiff, ydiff)
    # 6. Re-computes the entire 'nnbrs' matrix (an option slower than of 4. and 5. combined)
    
@jit   # FUNCTION 1
def find_diffusing_atom (rho: float, K_DIFF: np.array) -> (int, int):
    
    x0, y0 = 0, 0
    partial_sum = K_DIFF[x0, y0]
    
    while (rho > partial_sum):   # Continue adding (ie moving along 'K_DIFF') and checking
        
        if (x0 < LX - 1):   # Move along a 'K_DIFF' row
            x0 += 1
        else:   # The row has finished, move up
            x0 = 0
            y0 += 1
            
        partial_sum += K_DIFF[x0, y0]
    
    return x0, y0

@jit   # FUNCTION 2
def into_PBC (coord: int, L: int) -> int:
    
    coord = int(coord - L * np.floor(coord/L))
    
    return coord

@jit   # FUNCTION 3
def count_nbrs(x: int, y: int, height: np.array) -> int:
    
    n = 0
    z = height[x, y]
    
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
def update_diff_nnbrs (x0: int, y0: int, xdiff: int, ydiff: int, 
                       height: np.array, nnbrs: np.array) -> np.array:
    
    # Updates each position using the related deposition function.
    # It is somehow redundant because (x0, y0) and (xdiff, ydiff) are updated twice,
    # but I am lazy in working out a proper and slightly faster solution...
    
    nnbrs = update_depo_nnbrs(x0, y0, height, nnbrs)
    nnbrs = update_depo_nnbrs(xdiff, ydiff, height, nnbrs)
    
    return nnbrs

@jit   # FUNCTION 6
def update_nnbrs(height: np.array) -> np.array:
    
    nnbrs = np.zeros_like(height)
    
    for x in range(0, LX):
        for y in range(0, LY):
            nnbrs[x, y] = count_nbrs(x, y, height)
            
    return nnbrs

#%% -------------------- MAIN CODE -------------------- %%#

# ----- FUNDAMENTAL VARIABLES: initializations

f_output = read_input()

RNG = np.random.default_rng(seed=SEED)
height = np.zeros((LX, LY), dtype = np.int_)
sim_time = 0.   # In s
theta = 0.   # Actual measured coverage in ML/s or atom/s

# In the beginning all atoms are surface atoms (height == 0 everywhere),
# therefore each of them has 4 first neighbours.
#
# NOTE: 'height' can be < 0, but the algorithm still treats diffusion correctly,
#       that's because 'count_nbrs' operates also with negative numbers.

K_DIFF = np.ones((LX, LY), dtype = np.float_) * 4 * NU * np.exp( -(J0 + 4*J1)/(K_BOLTZ*TEMP) )

nnbrs = np.ones((LX, LY), dtype = np.int_) * 4

# ----- OUTPUT RELATED VARIABLES

snap_dt = (100 * 1e-3) * (0.2/PHI)   # In seconds
snap_time = 0.
n_snap = 0

# ----- MC LOOP: DEPOSITION ONLY ---> governed equivalently by 'sim_time' and 'theta'

while sim_time < END_TIME and DIFFUSION_BOOL == False:
    
    tau = - np.log(RNG.random()) / K_DEPO   # Escaping time in s, K_TOT is given only by K_DEPO
    
    # ----- DEPOSITION EVENT
    
    xdepo = RNG.integers(0, LX, endpoint=False)
    ydepo = RNG.integers(0, LY, endpoint=False)
    height[xdepo, ydepo] += 1
    
    # ----- END OF STEP ROUTINES: update 'sim_time' and compute/save outputs
    
    sim_time += tau
    
    # Regular output (at each step)
    
    theta = np.mean(height)   # Actual coverage
    rough = np.std(height)   # Actual roughness
    
    f_output.write(f"{sim_time:.9E} \t {PHI*sim_time:.9E} \t {tau:.9E} \t {theta:.9E} \t {rough:.9E} \n")
    
    # Save a snap of the atomic configuration every 'snap_dt'
    
    if (sim_time > snap_time):
        save_snap("SNAP" + str(n_snap).zfill(4) + ".xyz", theta, height)
        
        n_snap += 1
        snap_time += snap_dt
        
        print(f"{sim_time:E} / {END_TIME:E}")
        
# ----- MC LOOP: DEPOSITION AND DIFFUSION ---> governed by 'theta'

while theta < THETA and DIFFUSION_BOOL == True:
    
    summed_K_DIFF = np.sum(K_DIFF)   # Used in the output
    K_TOT = K_DEPO + summed_K_DIFF
    
    tau = - np.log(RNG.random()) / K_TOT   # Escape time in s
    
    # ----- CHOOSE THE EVENT
    
    rho = K_TOT * RNG.random()
    
    if (rho < K_DEPO):   # ----- DEPOSITION EVENT
        
        xdepo = RNG.integers(0, LX, endpoint=False)
        ydepo = RNG.integers(0, LY, endpoint=False)
        height[xdepo, ydepo] += 1
        
        # Update the 'nnbrs' matrix
        nnbrs = update_depo_nnbrs(xdepo, ydepo, height, nnbrs)
    
    else:   # ----- DIFFUSION EVENT
        
        rho -= K_DEPO
    
        x0, y0 = find_diffusing_atom(rho, K_DIFF)
        diff_dir = RNG.random()
        
        if (diff_dir < 0.25):   # East
            xdiff, ydiff = into_PBC(x0 + 1, LX), y0
            
        elif (diff_dir < 0.5):   # South
            xdiff, ydiff = x0, into_PBC(y0 - 1, LY)
            
        elif (diff_dir < 0.75):   # West
            xdiff, ydiff = into_PBC(x0 - 1, LX), y0
            
        else:   # North
            xdiff, ydiff = x0, into_PBC(y0 + 1, LY)
            
        height[x0, y0] -= 1
        height[xdiff, ydiff] += 1
        
        # Update the 'nnbrs' matrix
        nnbrs = update_diff_nnbrs(x0, y0, xdiff, ydiff, height, nnbrs)
    
    # ----- UPDATE THE 'K_DIFF' MATRIX
    
    E_b = J0 + nnbrs*J1
    K_DIFF = 4 * NU * np.exp( - E_b / (K_BOLTZ*TEMP) )
    
    # ----- END OF STEP ROUTINES: update 'sim_time' and compute/save outputs
    
    sim_time += tau
    
    # Regular output (at each step)
    
    theta = np.mean(height)   # Actual coverage
    rough = np.std(height)   # Actual roughness
    
    f_output.write(f"{sim_time:.9E} \t {PHI*sim_time:.9E} \t {tau:.9E} \t {theta:.9E} " \
                    f"\t {rough:.9E} \t {K_DEPO:.9E} \t {summed_K_DIFF:.9E} \n")
    
    # Save a snap of the atomic configuration every 'snap_dt'
    
    if (sim_time > snap_time):
        save_snap("SNAP" + str(n_snap).zfill(4) + ".xyz", theta, height)
        
        n_snap += 1
        snap_time += snap_dt
        
        print(f"{theta:f} / {THETA:f}")
    
# ----- END OF SIMULATION ROUTINES: things to do before ending code execution

f_output.close()

with open(OUTPUT_FOLDER + "/YOURInput_parameters.txt", "w") as input_params:
    input_params.write(str(DIFFUSION_BOOL) + "\n")
    input_params.write(str(PHI) + "\n")
    input_params.write(str(TEMP) + "\n")
    input_params.write(str(SEED) + "\n")
    input_params.write(str(snap_dt) + "\n")
    
#%% PROVE

f_output.close()