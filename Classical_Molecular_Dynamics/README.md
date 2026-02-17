For information about packages, start from here:

https://stackoverflow.com/questions/15746675/how-to-write-a-python-module-package

The code is intended to work only with Silver atoms interacting via the Lennard-Jones potential. In a future, and improbable, development, the functions 'LJ_6_12_energy', 'LJPOLY_energy' and similar SHOULD BE replaced with classes like 'Interaction_potential' and the relative methods, like 'Compute_energy' or even 'Compute_force'. The obtained structure should then be more suitable for different potential, and even different atoms if a class 'My_simulation' is added.

In short, a structure similar to the 'zpic' code must be envisioned.

## NOTE: launching several sims automatically
you can launch your simulations all at once using the 'launch_several_sims.py' script. It MUST be launched from the Anaconda Powershell Prompt, otherwise, inside the Spyder IDE, the std.output is not redirected toward the IPython console. If you want to use the Spyder IDE, you should activate the temporary status.txt file in 'main.py'.

## NOTE: changing machine
when you change the machine, i.e. the file system, you need to change the following:
- _**full_path**_ variable in the _**main.py**_ script
- _**full_path**_ variable in the _**my_module.py**_ script
- _**read_input**_ function in the _**my_module.py**_ script (not mandatory, see ***)
- _**folder_path**_ variable in the _**PostProcessing.py**_ script
- _**full_path**_ variable in the _**MinimizeConfiguration.py**_ script

Everything else will be managed by the program. 

# FILE SYSTEM STRUCTURE

The structure is defined as:

Folder_with_program
	|_____ my_package
	            |________ __pycache__
	            |________ "__init__.py"
	            |________ "my_module.py"
	            |________ "PostProcessing.py"
	            |________ "MinimizeConfiguration.py"
	            |________ "README.txt"
	|_____ Working_folder
	            |________ L10_22Ottobre
	            |________ L11_28Ottobre
	            |________      ...
			            |________ Generic_output_folder
							|________________ "Conf001.xyz"
							|________________ "Conf002.xyz"
							|________________ ...
							|________________ "Output.txt"
							|________________ "YOURInput_parameters.txt"
			            |________ "Input_parameters.txt"
			            |________ "launch_several_sims.py"
			            |________ "main.py"
	|_____ Useful_files_for_simulations
	            |________ "fcc100a108.txt"
	            |________ "fcc100a256.txt"
	            |________ ...
	            |________ "MINIMIZED_fcc111a336+1.txt"

*** The code was developed during course lectures, thus I had several lessons to deal with. 
    This useless nesting can be avoided and the 'lecture step' can be removed. Doing so, you
    must also change 'read_input', so that WORKING_DIR is not used anymore.

# INPUT file structure

The "Input_parameters.txt" file MUST list:

- name of the lattice file (string expected)
- three 0s or 1s deciding if and where PBCs must be used
- initial temperature (float expected)
- seed of the random number generator (int expected)
- boolean deciding whether to use the polynomial junction or not (use if == True)
- a temporary variable, often used for manually setting the time step size
- name of the output folder

An example is the following list:

fcc100a256.txt
1 1 0
20
16324478
True
2
Generic_output_folder
