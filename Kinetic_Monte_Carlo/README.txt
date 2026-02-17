The Kinetic Monte Carlo code is split in 2 different independent codes. No personal modules and packages are used, 
only those from the standard Python's library are present. The essential files are:

- 'main.py': responsibly for running the simulation.

- 'post_processing.py': analyze data and plot some graphs.

- '1st_part' and '2nd_part' folders, they were used in order to keep separated two different KMC exercises.
  If you want to get rid of them, just change the behaviour of 'read_input()' in the 'main.py' file.

- 'Input_parameters.txt': essential values to be passed for the simulation to work as defined by the user, otherwise
  the code itself will set the following values:
     - DIFFUSION_BOOL = False
     - PHI = 0.2
     - TEMP = 100.
     - SEED = 312428270054674116253981409899306192348
     - OUTPUT_FOLDER = "1st_part"

---------- The structure of the code is as follows:

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
					    |________________ "SNAP0000.xyz"
					    |________________ "SNAP0001.xyz"
					    |________________ ...
					    |________________ "Output.txt"
					    |________________ "YOURInput_parameters.txt"
	|________ "Input_parameters.txt"
	|________ "main.py"
	|________ "post_processing.py"

---------- The structure of the output from 'main.py' is clarified here below:

The folder storing the entirety of the simulation is called using a time convention. Example: '01_Jan31_10_35_phi_0.1'

'01': the month ordinal number (given in 01-12 format, used for file ordering)
'Jan': the abbreviated name of the month
'31': day of the month
'10': hour of the day (given in 00-23 format)
'35': minute of the hour (given in 00-59 format)
'0.1': the given deposition rate in units of [ML/s], that is the 'phi' cited in the filename

- REGULAR OUTPUT: escape time

The file is called 'Output.txt' and is divided in 5 columns:

0-th: 'sim_time', the length of the simulation until now. Given in seconds
1-th: the nominal coverage corresponding to 'sim_time', thus it is NOT the actual coverage. Given in MonoLayers (or atoms)
2-th: 'tau', the extracted escaping time of the selected (thus performed) event. Given in seconds
3-th: 'theta', the actual coverage. Given in MonoLayers (or atoms)
4-th: 'rough', the actual surface roughness. Given in MonoLayers (or atoms)

For simulations in the '2nd_part' only (that is diffusion only) the code outputs two more columns:
5-th: 'K_DEPO', the deposition transition rate
6-th: 'summed_K_DIFF', the total diffusion transition rate (sum along K_DIFF)

- EVERY NOW and THEN OUTPUT: snap

The files are called 'SNAP0001.xyz'. They are compatipable with Ovito and they list:

0-th row: total number of the atoms = atoms ON the surface + atoms OF the surface
1-th row: 'theta', the actual coverage. Given in MonoLayers (or atoms)
2-th and following rows: the position of each atom, given as 'x, y, z', in unit of atomic positions, which is zero-based 
     (thus '0, 34, 1' means 1st atom on the x axis, 35th on the y and 2nd on the z)

- FINAL OUTPUT: input_parameters

The file is called "YOURInput_parameters.txt" and shows different useful infos:

0-th row: DIFFUSION_BOOL - 'False' = Deposition only and save files in "1st_part", otherwise for 'True'
1-th row: PHI - Deposition flux expressed in ML/s
2-th row: TEMP - Temperature at which diffusion happens
3-th row: SEED - seed used for the random number generator

4-th row: snap_dt - the simulation time passing between two successive snaps (useful for Ovito's visualizations)

NOTE: the first 4 rows are those actually used in "Input_parameters.txt" (in "KMC/Ag epitaxy") when launching simulations