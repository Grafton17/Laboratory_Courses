The Metropolis Monte Carlo code is split in 3 different independent codes. No personal modules and packages are used, only those from the standard Python's library are present. The essential files are:

- _main.py_: responsibly for running the simulation.

- _launch_several_sims.py_: as hinted by its name, it eases the working pipeline. Just list all the simulations and the correspondent input values that you want to use and the sims will be executed by their own.

- _post_processing.py_: analyze data and plot some graphs.

- _1st_part_ and _2nd_part_ folders, they were used in order to keep separated two different MMC exercises. If you want to get rid of them, just change the behaviour of 'red_input()' in the 'main.py' file.

- _Input_parameters.txt_: essential values to be passed for the simulation to work as defined by the user, otherwise the code itself will set the following values:
     - TEMP_BOOL = False
     - TEMP = 100.
     - N_ADATOMS = 25
     - SEED = 312428270054674116253981409899306192348
     - OUTPUT_FOLDER = "1st_part"

# Code structure overview

Folder_with_program
	       |________ 1st_part
			                |________ Generic_output_folder
					                              |________________ "SNAP0000.xyz"
					                              |________________ "SNAP0001.xyz"
					                              |________________ ...
					                              |________________ "Output.txt"
					                              |________________ "YOURInput_parameters.txt"
	       |________ 2nd_part
			                |________ Generic_output_folder
					                              |________________ "Output.txt"
					                              |________________ "YOURInput_parameters.txt"
	       |________ _Input_parameters.txt_
	       |________ _launch_several_sims.py_
	       |________ _main.py_
	       |________ _post_processing.py_

## Managing the ouput of 'main.py'

The name of the output folder is something like "20_0_60_1", meaning:

- 20 adatoms are used
- temperature of 0 K
- the lattice size is equal to 60 atoms
- the used seed is the number 1

__REGULAR OUTPUT__: escape time

The file is called 'Output.txt' and is divided in 3 columns:

0-th: 'i', i-th step of the MC loop (starts from 1, bc the initialized configuration is stored as 0-th step)
1-th: 'old_Energy', the energy of the i-th configuration, whether it was accepted or rejected
2-th: 'is_accepted', 1 == the trial move was accepted, 0 == the trial move was rejected

__EVERY NOW and THEN OUTPUT__ (removed if TEMP_BOOl == True): snap

The files are called 'SNAP0001.xyz'. They are compatipable with Ovito and they list:

0-th row: total number of the atoms = atoms ON the surface + atoms OF the surface
1-th row: 'i', i-th step of the MC loop (starts from 1, bc the initialized configuration is stored as 0-th step)
2-th and following rows: the position of each atom, given as 'x, y, z', in unit of atomic positions, which is zero-based 
     (thus '0, 34, 1' means 1st atom on the x axis, 35th on the y and 2nd on the z)

- __FINAL OUTPUT__: input_parameters

The file is called "YOURinput_parameters.txt" and shows different useful infos:

0-th row: TEMP_BOOL - 'False' = No Metropolis rejection rule and save files in "1st_part", otherwise for 'True'
1-th row: TEMP - temperature at which the system sits
2-th row: N_ADATOMS - the number of adatoms on the surface
3-th row: SEED - seed used for the number random generator

NOTE: these rows are those actually used in "Input_parameters.txt" (in "MMC/Ag_surface_equilibrium") when launching simulations
      Except for the 'output_name'
