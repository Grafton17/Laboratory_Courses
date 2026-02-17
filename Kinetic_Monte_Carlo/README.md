The Kinetic Monte Carlo code is split in 2 different independent codes. No personal modules and packages are used, only those from the standard Python's library are present. The essential files are:

- __main.py__: responsibly for running the simulation.

- __post_processing.py__: analyze data and plot some graphs.

- __1st_part__ and __2nd_part__ folders, they were used in order to keep separated two different KMC exercises. If you want to get rid of them, just change the behaviour of 'read_input()' in the 'main.py' file.

- __Input_parameters.txt__: essential values to be passed for the simulation to work as defined by the user, otherwise the code itself will set the following values:
     - _**DIFFUSION_BOOL**_ = False
     - _**PHI**_ = 0.2
     - _**TEMP**_ = 100.
     - _**SEED**_ = 312428270054674116253981409899306192348
     - _**OUTPUT_FOLDER**_ = "1st_part"

# Code structure overview

> Folder_with_program
> 
>> __1st_part__
>>  
>>> Generic_output_folder
>>>  
>>>> "SNAP0000.xyz"
>>>> 
>>>> "SNAP0001.xyz"
>>>> 
>>>> ...
>>>> 
>>>> "Output.txt"
>>>> 
>>>> "YOURInput_parameters.txt"
>>>> 
>> __2nd_part__
>> 
>>> Generic_output_folder
>>>  
>>>> "SNAP0000.xyz"
>>>> 
>>>> "SNAP0001.xyz"
>>>> 
>>>> ...
>>> 
>>>> "Output.txt"
>>>> 
>>>> "YOURInput_parameters.txt"
>>>> 
>> __Input_parameters.txt__
>> 
>>  __main.py__
>> 
>> __post_processing.py__  

# Managing the ouput of 'main.py'

The folder storing the entirety of the simulation is called using a time convention. Example: _01_Jan31_10_35_phi_0.1_

- _**01**_: the month ordinal number (given in 01-12 format, used for file ordering)
- _**Jan**_: the abbreviated name of the month
- _**31**_: day of the month
- _**10**_: hour of the day (given in 00-23 format)
- _**35**_: minute of the hour (given in 00-59 format)
- _**0.1**_: the given deposition rate in units of \[ML/s\], that is the _**phi**_ cited in the filename 

## __REGULAR OUTPUT__

The file is called _**Output.txt**_ and is divided in 5 columns:

- 0-th: _**sim_time**_, the length of the simulation until now. Given in seconds
- 1-th: the nominal coverage corresponding to _**sim_time**_, thus it is NOT the actual coverage. Given in MonoLayers (or atoms)
- 2-th: _**tau**_, the extracted escaping time of the selected (thus performed) event. Given in seconds
- 3-th: _**theta**_, the actual coverage. Given in MonoLayers (or atoms)
- 4-th: _**rough**_, the actual surface roughness. Given in MonoLayers (or atoms)

For simulations in the '2nd_part' only (that is diffusion only) the code outputs two more columns:
- 5-th: _**K_DEPO**_, the deposition transition rate
- 6-th: _**summed_K_DIFF**_, the total diffusion transition rate (summed along K_DIFF)

## __EVERY NOW and THEN OUTPUT__ 

The files are called _**SNAP0001.xyz**_. They are compatipable with Ovito and they list:

- 0-th row: total number of the atoms = atoms ON the surface + atoms OF the surface
- 1-th row: _**theta**_, the actual coverage. Given in MonoLayers (or atoms)
- 2-th and following rows: the position of each atom, given as 'x, y, z', in unit of atomic positions, which is zero-based 
     (thus '0, 34, 1' means 1st atom on the x axis, 35th on the y and 2nd on the z)

## __FINAL OUTPUT__

The file is called _**YOURinput_parameters.txt**_ and shows different useful infos:

- 0-th row: _**DIFFUSION_BOOL**_ - 'False' = Deposition only and save files in "1st_part", otherwise for 'True'
- 1-th row: _**PHI**_ - Deposition flux expressed in ML/s
- 2-th row: _**TEMP**_ - Temperature at which diffusion happens
- 3-th row: _**SEED**_ - seed used for the number random generator
- 4-th row: _**snap_dt**_ - the simulation time passing between two successive snaps (useful for Ovito's visualizations)

NOTE: the first 4 rows are those actually used in "Input_parameters.txt" (in "KMC/Ag epitaxy") when launching simulations
