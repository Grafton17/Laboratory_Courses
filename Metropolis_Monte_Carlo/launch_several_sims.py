# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:06:45 2026

@author: SS882604
"""

import subprocess
import os
import sys

#%% EXERCISE 2.1

# sim_bool = [False, False, False, False, False, False, 
#             False, False, False, False, False, False, 
#             False, False, False, False, False, False, 
#             False]

# sim_temp = [0, 0, 0, 0, 0, 0, 
#             0, 0, 0, 0, 0, 0, 
#             0, 0, 0, 0, 0, 0, 
#             0]

# sim_N = [ 2,  3,  4,  5,  6,  7,  
#           8, 9 , 10, 11, 12, 13, 
#           14, 15, 16, 17, 18, 19, 
#           20]

# sim_seeds = [3, 3, 3, 3, 3, 3, 
#               3, 3, 3, 3, 3, 3, 
#               3, 3, 3, 3, 3, 3, 
#               3]   # Allowable: 1, 2, 3, 4, 5

# sim_names = ["2_0_60_3", "3_0_60_3", "4_0_60_3", "5_0_60_3", "6_0_60_3", "7_0_60_3", 
#               "8_0_60_3", "9_0_60_3", "10_0_60_3", "11_0_60_3", "12_0_60_3", "13_0_60_3", 
#               "14_0_60_3", "15_0_60_3", "16_0_60_3", "17_0_60_3", "18_0_60_3", "19_0_60_3", 
#               "20_0_60_3"]

# # This list MUST NOT be changed, you modify the seed via 'sim_seeds', which is one-base

# seeds = [312428270054674116253981409899306192348,
#           144193831631783441215101441714363001219,
#           264914786653715399394480469990166609730,
#           332496105924005586042363788695332076602,
#           42179434837711941487287880505023012954]

#%% EXERCISE 2.2

# sim_bool = [True, True, True, True, True, True, 
#             True, True, True, True, True, True,
#             True]

# sim_temp = [200, 300, 600, 750, 900, 1100,
#             1200, 1300, 1400, 1600, 1750, 1900, 
#             2250]

# sim_N = [25, 25, 25, 25, 25, 25,  
#           25, 25, 25, 25, 25, 25,
#           25]

# sim_seeds = [1, 1, 1, 1, 1, 1,
#               1, 1, 1, 1, 1, 1,
#               1]   # Allowable: 1, 2, 3, 4, 5

# sim_names = ["25_200_60_1", "25_300_60_1", "25_600_60_1", "25_750_60_1", "25_900_60_1", "25_1100_60_1", 
#               "25_1200_60_1", "25_1300_60_1", "25_1400_60_1", "25_1600_60_1", "25_1750_60_1", "25_1900_60_1",
#               "25_2250_60_1"]

# sim_bool = [True, True, True, True, True, True, 
#             True, True, True, True, True, True, 
#             True, True, True, True, True, True, 
#             True, True, True, True]

# sim_N = [25, 25, 25, 25, 25, 25,  
#           25, 25, 25, 25, 25, 25,
#           25, 25, 25, 25, 25, 25,
#           25, 25, 25, 25]

# sim_names = ["25_1_20_3", "25_10_20_3", "25_100_20_3", "25_200_20_3", "25_300_20_3", "25_500_20_3", 
#               "25_600_20_3", "25_750_20_3", "25_900_20_3", "25_1000_20_3", "25_1100_20_3", "25_1200_20_3",
#               "25_1300_20_3", "25_1400_20_3", "25_1500_20_3", "25_1600_20_3", "25_1750_20_3", "25_1900_20_3",
#               "25_2000_20_3", "25_2250_20_3", "25_2500_20_3", "25_2700_20_3"]

# sim_temp = [1, 10, 100, 200, 300, 500, 
#             600, 750, 900, 1000, 1100, 1200, 
#             1300, 1400, 1500, 1600, 1750, 1900, 
#             2000, 2250, 2500, 2700]

# sim_seeds = [3, 3, 3, 3, 3, 3,
#              3, 3, 3, 3, 3, 3,
#              3, 3, 3, 3, 3, 3,
#              3, 3, 3, 3]

# This list MUST NOT be changed, you modify the seed via 'sim_seeds', which is one-base

seeds = [312428270054674116253981409899306192348,
          144193831631783441215101441714363001219,
          264914786653715399394480469990166609730,
          332496105924005586042363788695332076602,
          42179434837711941487287880505023012954]

#%% SINGLE SIM

sim_bool = [True]
sim_temp = [1250]
sim_N = [25]
sim_seeds = [3]
sim_names = ["25_1250_40_3"]

#%% MAIN

input_file = r"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/MMC/Ag_surface_equilibrium/Input_parameters.txt"

for i in range(0, len(sim_names)):   # len(sim_names) 
    
    print("\nStarting", sim_names[i], end="\n\n")
    
    # Change 'Input_parameters.txt'
    f_input = open(input_file, 'w')
    
    input_text = str(sim_bool[i]) + f"\n{sim_temp[i]}\n{sim_N[i]}\n" + str(seeds[sim_seeds[i]-1]) + "\n" + sim_names[i]
    f_input.write(input_text)
    
    f_input.close()
    
    # Launch the simulation
    os.system("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/MMC/Ag_surface_equilibrium//main.py\"")
    
    print("############### Finished", sim_names[i], end=" ###############\n")
    
