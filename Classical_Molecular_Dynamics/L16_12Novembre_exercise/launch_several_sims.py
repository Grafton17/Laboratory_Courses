# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 22:12:28 2025

@author: SS882604
"""

import subprocess
import os
import sys

#%% EXERCISE 1

sim_names = [ "POLY_2000K_1fs_5000iter", "POLY_2000K_2fs_5000iter",
              "POLY_2000K_3fs_5000iter", "POLY_2000K_4fs_5000iter",
              "POLY_2000K_5fs_5000iter", "POLY_2000K_6fs_5000iter"
              ]

sim_dt = [1, 2, 3, 4, 5, 6]

sim_junction = [True, True, True, True, True, True]

Temp = int(2000)

#%% TESTS

# sim_names = ["test_1"]
# sim_dt = [15]
# sim_junction = [True]

#%% MAIN

input_file = r"Y:\My Drive\MIA_UNI\Magistrale\CMS\Working_folder\L16_12Novembre_exercise\Input_parameters.txt"

for i in range(0, len(sim_names)):   # len(sim_names) 
    
    print("\nStarting", sim_names[i], end="\n\n")
    
    # Change 'Input_parameters.txt'
    f_input = open(input_file, 'w')
    
    input_text = f"fcc100a256.txt\n{Temp}\n33282739\n{sim_junction[i]}\n{sim_dt[i]}\n{sim_names[i]}"
    f_input.write(input_text)
    
    f_input.close()
    
    # Launch the simulation
    os.system("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/main.py\"")
    
    print("############### Finished", sim_names[i], end=" ###############\n")
    

#%%
# subprocess.run("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/test.py\"", capture_output=True)
# os.system("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/test.py\"")

# subprocess.call("python \"Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/L16_12Novembre_exercise/test.py\"", 
#                stderr=sys.stderr, stdout=sys.stdout)