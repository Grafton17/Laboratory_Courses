# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 22:12:28 2025

@author: SS882604
"""

import subprocess
import os
import sys

#%% EXERCISE 1

sim_names = ["1100K_2fs_10ps", "1200K_2fs_10ps", "1300K_2fs_10ps", "1400K_2fs_10ps", "1500K_2fs_10ps",
             "1600K_2fs_10ps", "1700K_2fs_10ps", "1800K_2fs_10ps", "1900K_2fs_10ps"]

sim_dt = 2

sim_junction = True

Temp = [1300]

#%% TESTS

# sim_names = ["test_1"]
# sim_dt = [15]
# sim_junction = [True]

#%% MAIN

input_file = r"Y:\My Drive\MIA_UNI\Magistrale\CMS\Working_folder\L20-22_03Dicembre_exercise\Input_parameters.txt"

for i in range(0, len(sim_names)):   # len(sim_names) 
    
    print("\nStarting", sim_names[i], end="\n\n")
    
    # Change 'Input_parameters.txt'
    f_input = open(input_file, 'w')
    
    input_text = f"MINIMIZED_fcc100a256.txt\n0 0 0\n{Temp[i]}\n33282739\n{sim_junction}\n{sim_dt}\n{sim_names[i]}"
    f_input.write(input_text)
    
    f_input.close()
    
    # Launch the simulation
    os.system("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L20-22_03Dicembre_exercise/main.py\"")
    
    print("############### Finished", sim_names[i], end=" ###############\n")
    

#%%
# subprocess.run("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/test.py\"", capture_output=True)
# os.system("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/test.py\"")

# subprocess.call("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/test.py\"", 
#                stderr=sys.stderr, stdout=sys.stdout)
