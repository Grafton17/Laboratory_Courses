# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 22:10:36 2025

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/"   # Where the module sits
# folder names: "Working_folder/"   "Useful_files_for_simulations/"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.constants as cst
from numba import jit

import os
import sys

mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['figure.figsize'] = (8, 6)
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.bbox'] = 'tight'

# Global variables - MUST be in uppercase

MASS = 108 * 1.66e-27 / 16   # Silver atomic mass = 108 u.m.a.
EPS, SIGMA = 0.345, 2.644   # LJ_6_12 values for Silver in eV and angstrom
K_BOLTZ = 1/11603   # 1 eV per 11603 K of temperature

WORKING_DIR = ""   # Stores the name of the working directory
LATTICE_FILE = ""   # Stores the lattice's file's name

LJPOLY = ""   # Stores the polynomial junction
A, B, C, D, E, F, G, H = 0., 0., 0., 0., 0., 0., 0., 0.
RC = 4.5   # Cut-off radius in angstrom
RP = .0   # Polynomial junction radius in angstrom

L_PBC = np.array([0., 0., 0.])   # PBCs cell size along x, y and z
# L_LATTICE = np.array([16.6416, 16.6416, 16.6416])   # Lattice size for point 4
L_LATTICE = np.array([20.3817, 20.3817, 20.3817])   # Lattice size for point 5

TEMP = .0   # Temporary variable for several different things

# NOTE: if RP > RC, the polynomial is initialized but 'compute_Epot' and 
#       'compute_force' never use it, because they use only the LJ potential.
#       Thus, setting RP > RC essentially remove the presence of the junction.
#       This is done by passing 'False' as 'use_junction' in 'Input_parameters.txt'

# NOTE: if L_PBC[0] > 16.6416 by several times RC, then the PBCs are killed 
#       along the 0-th direction. The directions are ordered as x, y and z.
#       That's because the nearest image of the k-th atom to the j-th atom
#       would be the k-th atom itself.

#%% --------------- Class definitions ---------------

class poly7:
    def __init__(self, eps, sig, rcut, rpoly):
        
        self.A = ( (1/((rcut - rpoly)**7* rpoly**12))*4* eps* rcut**4* sig**6*
                 (2*rpoly**6 *(-42* rcut**3 + 182*rcut**2* rpoly -273* rcut* rpoly**2 +
          143* rpoly**3) + (455* rcut**3 - 1729*rcut**2* rpoly +
          2223* rcut* rpoly**2 - 969* rpoly**3)* sig**6) )
    
        self.B = ( (1/((rcut - rpoly)**7* rpoly**13))*16* eps* rcut**3* sig**6* 
                    (rpoly**6* (54* rcut**4 - 154* rcut**3* rpoly + 
          351* rcut *rpoly**3 - 286* rpoly**4) + (-315* rcut**4 + 749 *rcut**3 * rpoly + 
          171 *rcut**2* rpoly**2 - 1539* rcut* rpoly**3 + 969* rpoly**4)* sig**6) )
                                                
        self.C = ( (1/((rcut - rpoly)**7* rpoly**14))*12* eps* rcut**2* sig**6* 
          (rpoly**6* (-63* rcut**5 - 7* rcut**4 *rpoly + 665* rcut**3 *rpoly**2 - 
          975* rcut**2* rpoly**3 - 52* rcut* rpoly**4 + 572* rpoly**5) + 
          2 *(195* rcut**5 + 91* rcut**4* rpoly - 1781* rcut**3* rpoly**2 + 
          1995 *rcut**2* rpoly**3 + 399* rcut* rpoly**4 - 969* rpoly**5)* sig**6) )
              
        self.D = ( (1/((rcut - rpoly)**7* rpoly**15))*16* eps* sig**6* 
          (rcut* rpoly**6* (14* rcut**6 + 126* rcut**5* rpoly - 420* rcut**4* rpoly**2 
          -90* rcut**3* rpoly**3 + 1105* rcut**2* rpoly**4 - 624* rcut* rpoly**5 - 
          286 *rpoly**6) + rcut* (-91* rcut**6 - 819* rcut**5* rpoly + 2145 * 
          rcut**4 * rpoly**2 + 1125* rcut**3* rpoly**3 - 5035* rcut**2* rpoly**4 + 
          1881* rcut* rpoly**5 + 969* rpoly**6)* sig**6) )
                                 
        self.E = ( (1/((rcut - rpoly)**7* rpoly**15))*4* eps* sig**6* 
          (2* rpoly**6* (-112* rcut**6 - 63* rcut**5* rpoly + 1305* rcut**4* rpoly**2 
          -1625* rcut**3* rpoly**3 - 585* rcut**2* rpoly**4 + 
          1287 *rcut* rpoly**5 + 143* rpoly**6) + (1456*rcut**6 +1404*rcut**5* rpoly - 
          14580 *rcut**4* rpoly**2 + 13015* rcut**3* rpoly**3 + 7695* rcut**2* rpoly**4 - 
          8721 *rcut* rpoly**5 - 969* rpoly**6)* sig**6) )
                                                 
        self.F = ( (1/((rcut - rpoly)**7* rpoly**15))*48* eps* sig**6* 
          (-rpoly**6* (-28* rcut**5 + 63* rcut**4* rpoly + 65* rcut**3* rpoly**2 - 
          247* rcut**2* rpoly**3 + 117* rcut* rpoly**4 + 65* rpoly**5) + 
          (-182* rcut**5 + 312* rcut**4* rpoly + 475* rcut**3* rpoly**2 - 
          1140* rcut**2* rpoly**3 + 342* rcut* rpoly**4 + 228* rpoly**5)* sig**6) )
              
        self.G = ( (1/((rcut - rpoly)**7* rpoly**15))*4* eps* sig**6* (rpoly**6* 
          (-224* rcut**4 + 819* rcut**3* rpoly - 741* rcut**2* rpoly**2 
          -429* rcut* rpoly**3 + 715* rpoly**4) + 2 *(728* rcut**4 
          -2223* rcut**3* rpoly + 1425* rcut**2* rpoly**2 
          +1292* rcut* rpoly**3 - 1292* rpoly**4)* sig**6) )
                              
        self.H = ( (1/((rcut - rpoly)**7* rpoly**15))*16* eps* sig**6* (rpoly**6* 
          (14* rcut**3 - 63* rcut**2* rpoly + 99* rcut* rpoly**2 - 55* rpoly**3) + 
          (-91* rcut**3 + 351* rcut**2* rpoly - 459* rcut* rpoly**2 
          +204* rpoly**3)* sig**6) )
        
    def evaluate(self, r):
        
        evaluation = (self.A        + self.B * r    + 
                      self.C * r**2 + self.D * r**3 + 
                      self.E * r**4 + self.F * r**5 +
                      self.G * r**6 + self.H * r**7)
        
        return evaluation
    
    def evaluate_1Der(self, r):
        
        evaluation = (self.B            + self.C * 2 * r    + 
                      self.D * 3 * r**2 + self.E * 4 * r**3 +
                      self.F * 5 * r**4 + self.G * 6 * r**5 +
                      self.H * 7 * r**6)
        
        return evaluation

#%% --------------- I/O Functions ---------------

# The 'Input_parameters.txt' file must be composed as:
    # 'LATTICE_FILE' name
    # 'use_PBCs': a string made of three 0s or 1s separated by white spaces.
    #             For example "1 0 1", kills PBCs along y and not along x and z
    # 'Temp' value in Kelvin
    # 'myseed' value used for random extractions, must be an INTEGER
    # 'use_junction' a bool deciding whether to use the polynomial junction or not
    # 'TEMPORARY' used differently in different simulations
    # 'output_folder' name (where you want the results of the simulation for post-processing)
    
# Output values:
    # 'pos': matrix of the atoms' coords in angstrom
    #        pos[i, 2] indicates the z coord of the i-th atom
    # 'Temp' same as before
    # 'myseed' an int number
    # 'output_folder' a string

def read_input (input_file = "Input_parameters.txt"):
    
    # ----- Find the working directory (ie the last lesson) ------
        # This mess could have been avoided if I had noticed the 
        # 'Browse working directory' button in the upper right corner of Spyder...
        # Nevermind, I'll remember it for the next time -_-
    
    global WORKING_DIR
    
    # dirs = os.listdir(full_path + "Working_folder")
    # dirs.sort()
    # WORKING_DIR = dirs[-1]   # Always refer to the last lesson in the "Working_folder"
    
    WORKING_DIR = "L20-22_03Dicembre_exercise"
    
    # ----- Read the input parameters from the 'input_file' ------
    
    filename = full_path + "Working_folder/" + WORKING_DIR + "/" + input_file
    
    with open(filename, 'r') as input_params:
    
        global LATTICE_FILE
        global L_PBC
        global L_LATTICE
        global RP
        global RC
        global TEMPORARY
    
        LATTICE_FILE = input_params.readline()
        LATTICE_FILE = LATTICE_FILE.splitlines()[0]   # Remove the '\n' at the end
        
        use_PBCs = input_params.readline()
        use_PBCs = use_PBCs.split(maxsplit = 2)
        
        L_PBC[0] = L_LATTICE[0] + (1 - int(use_PBCs[0])) * 20*RC
        L_PBC[1] = L_LATTICE[1] + (1 - int(use_PBCs[1])) * 20*RC   # Note: (1 - int(...)) switches 0 and 1
        L_PBC[2] = L_LATTICE[2] + (1 - int(use_PBCs[2])) * 20*RC
        
        Temp = float(input_params.readline())
        myseed = int(input_params.readline())
        
        use_junction = input_params.readline()
        use_junction = use_junction.splitlines()[0]
        
        if use_junction == "True":
            RP = 4.2
        else:
            RP = RC + 0.3
            
        TEMPORARY = float(input_params.readline())
        
        output_folder = input_params.readline()
        output_folder = output_folder.splitlines()[0]
        output_folder = full_path + "Working_folder/" + WORKING_DIR + "/" + output_folder   # Absolute path needed
        
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            if (input("The output folder already exists. Rewrite it [y/N]?\t") == "N"):
                print("\n")
                raise 
            
    global LJPOLY
    global A, B, C, D, E, F, G, H
            
    LJPOLY = poly7(EPS, SIGMA, RC, RP)
    
    A, B, C, D = LJPOLY.A, LJPOLY.B, LJPOLY.C, LJPOLY.D
    E, F, G, H = LJPOLY.E, LJPOLY.F, LJPOLY.G, LJPOLY.H
    
    # ----- Upload the lattice -----
    
    pos_file = full_path + "Useful_files_for_simulations/" + LATTICE_FILE
    pos = np.loadtxt(pos_file, skiprows=2)   # Atoms' positions
    
    # print(f"Uploaded lattice: '{LATTICE_FILE}'\n")
    
    return pos, Temp, myseed, output_folder

#%% --------------- Lennard-Jones 6-12 pair potential ---------------

@jit
def LJ_6_12_energy (dist: float) -> float:
    
    # NOTE: 'dist' is assumed to be given in angstroms
    # NOTE: 'EPS' is assumed to be given in eV
    # NOTE: 'SIGMA' is asssumed to be given in angstroms
    
    frac = SIGMA/dist
    Epot = 4*EPS*(frac**12 - frac**6)
    
    return Epot

@jit
def LJ_6_12_force (dist: float, x_fix: np.ndarray, x_ext: np.ndarray) -> np.ndarray:
    
    # NOTE: 'dist' is assumed to be given in angstroms
    # NOTE: 'EPS' is assumed to be given in eV
    # NOTE: 'SIGMA' is asssumed to be given in angstroms
    
    # The force due to atom 'ext' acts on atom 'fix'
    
    delta = x_fix - x_ext
    delta = delta - L_PBC * np.round(delta/L_PBC)
    S = SIGMA**6
    force = 24*EPS*S * (1/dist)**8 * delta * (2*S * (1/dist)**6 - 1)
    
    return force

@jit
def LJPOLY_energy (dist: float) -> float:
    
    # NOTE: 'dist' is assumed to be given in angstroms
    
    Epot = (A           + B * dist    + 
            C * dist**2 + D * dist**3 + 
            E * dist**4 + F * dist**5 +
            G * dist**6 + H * dist**7)
    
    return Epot

@jit
def LJPOLY_force (dist: float, x_fix: np.ndarray, x_ext: np.ndarray) -> np.ndarray:
    
    # NOTE: 'dist' is assumed to be given in angstroms
    
    # The force due to atom 'ext' acts on atom 'fix'
    
    delta = x_fix - x_ext
    delta = delta - L_PBC * np.round(delta/L_PBC)
    der_LJPOLY = (B               + C * 2 * dist    + 
                  D * 3 * dist**2 + E * 4 * dist**3 + 
                  F * 5 * dist**4 + G * 6 * dist**5 +
                  H * 7 * dist**6)
    
    force = - der_LJPOLY * delta / dist
    
    return force

#%% --------------- Verlet list, force and potential energy ---------------

# NOTE: 'nnbrs', 'whichnbr', 'force' and 'Epot' must be initialized to zero

@jit
def compute_Verlet (whichnbr: np.ndarray, pos: np.ndarray, natoms: int) -> list[np.ndarray]:
    
    nnbrs = np.zeros(natoms, dtype=np.int_)   # np.int_ is needed for @jit
    
    for i in range(0, natoms):
        for j in range(0, natoms):
            
            if i == j:
                continue
            
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            dx = dx - L_PBC[0] * np.round(dx/L_PBC[0])
            dy = dy - L_PBC[1] * np.round(dy/L_PBC[1])
            dz = dz - L_PBC[2] * np.round(dz/L_PBC[2])
            dist2 = dx*dx + dy*dy + dz*dz
            dist = dist2**0.5
            
            if dist > RC:   # When dist > RC, the potential is null
                continue
            
            # nnbrs[i] starts at 0 and gets incremented as neighbours are found
            
            whichnbr[i, nnbrs[i]] = j
            nnbrs[i] += 1

    # print("Done - VERLET\n")
    
    return nnbrs, whichnbr

@jit
def compute_fastVerlet (whichnbr: np.ndarray, pos: np.ndarray, natoms: int) -> list[np.ndarray]:
    
    nnbrs = np.zeros(natoms, dtype=np.int_)   # np.int_ is needed for @jit
    
    for i in range(0, natoms - 1):
        for j in range(i + 1, natoms):
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            dx = dx - L_PBC[0] * np.round(dx/L_PBC[0])
            dy = dy - L_PBC[1] * np.round(dy/L_PBC[1])
            dz = dz - L_PBC[2] * np.round(dz/L_PBC[2])
            dist2 = dx*dx + dy*dy + dz*dz
            dist = dist2**0.5
            
            if dist > RC:   # When dist > RC, the potential is null
                continue
            
            whichnbr[i, nnbrs[i]] = j
            nnbrs[i] += 1
            
            whichnbr[j, nnbrs[j]] = i
            nnbrs[j] += 1
            
    # print("Done - fastVERLET\n")
    
    return nnbrs, whichnbr

    # FOR CHECKING the two VERLET functions, I USED:
    # a = np.equal(whichnbr, whichnbr_fast)   # Confront two arrays element wise
    # are_they_entirely_equal = np.all(a)   # Test whether all a's elements are 'True' or not

@jit
def compute_force (nnbrs: np.ndarray, whichnbr: np.ndarray, pos: np.ndarray, natoms: int) -> np.ndarray:
    
    force = np.zeros((natoms, 3))
    
    for i in range(0, natoms):
        for j in range(0, nnbrs[i]):
            k = whichnbr[i, j]
            
            dx = pos[i, 0] - pos[k, 0]
            dy = pos[i, 1] - pos[k, 1]
            dz = pos[i, 2] - pos[k, 2]
            dx = dx - L_PBC[0] * np.round(dx/L_PBC[0])
            dy = dy - L_PBC[1] * np.round(dy/L_PBC[1])
            dz = dz - L_PBC[2] * np.round(dz/L_PBC[2])
            dist2 = dx*dx + dy*dy + dz*dz
            dist = dist2**0.5
            
            if (dist < RP):   # Before the polynomial junction
                force[i] += LJ_6_12_force(dist, pos[i], pos[k])
            else:
                force[i] += LJPOLY_force(dist, pos[i], pos[k])

    # print("Done - FORCE\n")
    
    return force

@jit
def compute_Epot (nnbrs: np.ndarray, whichnbr: np.ndarray, pos: np.ndarray, natoms: int) -> float:
    
    Epot = 0    
    
    for i in range(0, natoms):
        for j in range(0, nnbrs[i]):
            k = whichnbr[i, j]
            
            dx = pos[i, 0] - pos[k, 0]
            dy = pos[i, 1] - pos[k, 1]
            dz = pos[i, 2] - pos[k, 2]
            dx = dx - L_PBC[0] * np.round(dx/L_PBC[0])
            dy = dy - L_PBC[1] * np.round(dy/L_PBC[1])
            dz = dz - L_PBC[2] * np.round(dz/L_PBC[2])
            dist2 = dx*dx + dy*dy + dz*dz
            dist = dist2**0.5
            
            if (dist < RP):   # Before the polynomial junction
                Epar = LJ_6_12_energy(dist)   # Partial energy
            else:
                Epar = LJPOLY_energy(dist)
            
            Epot += Epar

    Epot = Epot/2
    
    # print("Done - ENERGY\n")
    
    return Epot

#%% --------------- Kinetic energy and velocity initialization ---------------

# In this MD code, we extract initial velocities from a uniform distribution (PDF)
# The boundaries of PDF, which is 0 centered, are computed in order to have:
#      1/2 * m * <v**2> =  1/2 * k_B * T, where <v**2> is a PDF average

# The number of extractions is not infinite, thus it generally happens that:
#  - <v> != 0
#  - <v**2> != k_B * T / m
# Therefore, for each direction, we shift the PDF by <v> and renormalize it
# so that <v> = 0 (the lattice does not translate) and that <v**2> (T is 
# correctly restored in the average sense of the micro-canonical ensemble).

# Moreover, the renormalization accounts also for the effects of PDF shifting
# on the value of <v**2>.

# Read 'TsdM5En.pdf' for more.

def compute_Ekin (vel: np.ndarray, natoms: int) -> float:
    
    Ekin = 0
    
    for i in range(0, natoms):
        for j in range(0,3):
            Epar = 0.5 * MASS * vel[i,j]**2   # Partial energy
            Ekin += Epar
    
    # print("Done - KINETIC\n")
    
    return Ekin

def compute_fastEkin (vel: np.ndarray) -> float:
    
    Ekin = 0.5 * MASS * np.sum(vel**2)
    
    return Ekin

def compute_initialVelocity (Temp: float, natoms: int, myseed: int) -> np.ndarray:
    
    # Initializations
    vel = np.zeros((natoms, 3))   # Matrix of the velocities
    np.random.seed(myseed)   # Set the seed for the random number generator
    CC = np.sqrt(3*K_BOLTZ*Temp/MASS)   # Boundaries of the uniform distribution
    
    for i in range(0, natoms):
        for j in range(0,3):
            vel[i, j] = np.random.uniform(-CC, CC)
            
    # Fix the fact that <v> != 0 for a finite lattice
    vel[:, 0] = vel[:, 0] - np.mean(vel[:,0])
    vel[:, 1] = vel[:, 1] - np.mean(vel[:,1])
    vel[:, 2] = vel[:, 2] - np.mean(vel[:,2])
    
    # Fix the fact that <v**2> is not what is expected from an infinite lattice
    Ekin = compute_fastEkin(vel)
    Temp_extracted = Ekin / (natoms * K_BOLTZ * 1.5)
    vel = vel * np.sqrt(Temp/Temp_extracted)
    
    # print("Done - initialVELOCITY\n")
    
    return vel

#%% --------------- FANCY DRAWINGS ---------------

# Plot the lattice with an atom and its nbrs highlighted

def plot_nbrs (atom, whichnbr, nnbrs, pos, natoms):
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(pos[:,0], pos[:,1], pos[:,2], marker='o', color='gray', s = 300)

    for i in range(0, nnbrs[atom]):
        k = whichnbr[atom, i]
        ax.scatter(pos[k,0], pos[k,1], pos[k,2], marker='o', color='k', s = 300)
    
    ax.scatter(pos[atom, 0], pos[atom, 1], pos[atom, 2], marker='o', color='r', s = 300)
        
    ax.set_xlabel('X coord')
    ax.set_ylabel('Y coord')
    ax.set_zlabel('Z coord')

    plt.show()
    
# Plot the velocity distribution for the x, y and z directions
    
def plot_velDistribution (vel, bins=10):
    
    plt.hist(vel[:, 0], alpha=0.5, label="x", bins=bins)
    plt.hist(vel[:, 1], alpha=0.5, label="y", bins=bins)
    plt.hist(vel[:, 2], alpha=0.5, label="z", bins=bins)
    plt.legend()
    plt.grid()
    plt.title("Velocity distribution")
    plt.ylabel("Absolute counts")
    
    #Check: https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#:~:text=Boltzmann%20speed%20distribution-,Relaxation%20to%20the%202D%20Maxwell%E2%80%93Boltzmann%20distribution,-%5Bedit%5D

    plt.show()
    
# Plot the (expected) Maxwell-Boltzmann distribution
    
def plot_speedDistribution (vel, bins=10):
    
    # Some useful values and manipulations
    natoms = np.shape(vel)[0]
    
    vel_sq = np.sum(vel**2, axis=1) 
    speed = np.sqrt(vel_sq)
    
    Temp = 0.5 * MASS * np.sum(vel_sq) / (natoms * K_BOLTZ * 1.5)

    # Maxwell-Boltzmann expected distribution
    MB = lambda vv: (natoms * (MASS * 16 / (2 * cst.pi * cst.k * Temp) )**1.5 *
                    4 * cst.pi * (vv*1e-10)**2 * 
                    np.exp(- MASS * vv**2 / (2 * K_BOLTZ * Temp) ) )
    
    # Plotting stuff
    max_v = np.max(speed)
    x_axis = np.linspace(0, max_v, 1000)

    plt.hist(speed, bins, label='Sampled')
    plt.plot(x_axis, MB(x_axis), label='Expected')
    
    plt.grid()
    plt.title("Maxwell-Boltzmann distribution")
    plt.xlabel("Speed")
    plt.ylabel("Absolute counts")
    plt.legend()
    
    plt.show()

#%% --------------- POST-PROCESSING ---------------

# Plot the time evolution of the total energy and of the instantaneous temperature

def plot_timeEvolution (time: np.ndarray, evolution_Etot: np.ndarray, 
                        evolution_T: np.ndarray, useOffset: bool = False) -> None:
    """Plot ('total energy' vs 'time') and ('temperature' vs 'time')."""
    
    fig, [axE, axT] = plt.subplots(2, 1, figsize=(8,6), layout = 'constrained')

    axE.plot(time, evolution_Etot)
    axE.ticklabel_format(useOffset=useOffset)
    axE.grid()
    axE.set_title("Total energy")
    axE.set_ylabel("Energy [eV]")

    axT.plot(time, evolution_T)
    axT.grid()
    axT.set_title("Temperature")
    axT.set_ylabel("Temperature [K]")
    
    fig.supxlabel("time [s]")
    
    plt.show()
    
    return fig, axE, axT
    
def compute_RelFluct (data: np.ndarray, time: np.ndarray, thermal_time: float = 3 * 1e-12, 
                      verbose = True) -> (float, float):
    """Compute the mean, absolute and relative fluctuation of 'data' after thermalization has occured."""
    
    thermal_index = 0

    for i in range(0, len(data)):
        if (time[i] > thermal_time ):
            thermal_index = i
            break
        
    data_mean = np.mean(data[thermal_index:])

    data2_mean = np.mean(data[thermal_index:]**2)
    fluct = np.sqrt(data2_mean - data_mean**2)
    # E_fluct = np.std(evolution_Etot, mean=E_mean)   # <--- Available in newer Numpy versions

    rel_fluct = fluct / np.abs(data_mean)
    
    if verbose:
        print(f"<S>: {data_mean} \t Std: {fluct} \t rel_fluct: {rel_fluct:.3E}")
    
    return data_mean, fluct, rel_fluct

def compute_InCellTraj (traj: np.ndarray, pbc_size: list[float]) -> np.ndarray:
    """Compute the trajectory of the adatom using the PBC."""
    
    pbc_size = np.array(pbc_size)
    cell_traj = traj - pbc_size * np.floor(traj/pbc_size)
    
    return cell_traj

def toPBC_ConfFile (old_folder: str = "point4_1700K_2fs_300ps") -> None:
    """Convert the 'Conf###.xyz' file into a new one where the adatom trajectory is given in PBCs."""
    
    # Obtain the name of the new_conf folder and eventually create it
    splitted = old_folder.rsplit(sep = "_", maxsplit = 3)
    new_folder = splitted[0] + "_" + splitted[1] + "_PBC"
    new_path = full_path + "Working_folder/L20-22_03Dicembre_exercise/" + new_folder
    
    try:
        os.mkdir(new_path)
    except FileExistsError:
        if (input("'" + new_folder + "' folder already exists. Rewrite it [y/N]?\t") == "N"):
            print("\n")
            raise
            
    # Scan the conf files, convert them into PBC trajectories and save them in the new_conf folder
    old_path = full_path + "Working_folder/L20-22_03Dicembre_exercise/" + old_folder
    
    n_Conf = len(os.listdir(old_path)) - 2   # Number of 'Conf###.xyz' files
    
    if (splitted[0] == "point4"):   # fcc100a256+1 lattice
        pbc_size = [16.6416, 16.6416, 1000]
    elif (splitted[0] == "point5"):   # fcc111a336+1 lattice
        pbc_size = [20.3817, 20.3817, 1000]
    else:
        print("Not found the correspondent 'pbc_size' list\n")
        sys.exit()
    
    for i in range(1, n_Conf + 1):
        old_file = old_path + "/Conf" + str(i).zfill(3) + ".xyz"
        data = np.loadtxt(old_file, skiprows = 2)
        data[-1,:] = compute_InCellTraj(data[-1,:], pbc_size)   # Last row: adatom's coords
        
        new_file = new_path + "/Conf" + str(i).zfill(3) + ".xyz"
        np.savetxt(new_file, data, header = f"{np.size(data, axis=0)}\nTime={i*200*1e-15}", comments="")
        
        if (i % 50 == 0):
            print(i)
        
    return None
       
def extract_AdatomTraj (coords_folder: str = "point4_1700K_PBC") -> str:
    
    coords_path = full_path + "Working_folder/L20-22_03Dicembre_exercise/" + coords_folder
    
    n_Conf = len(os.listdir(coords_path)) - 2   # Number of 'Conf###.xyz' files
    ex_name = coords_folder.rsplit(sep = "_", maxsplit = 3)[0]
    
    if (ex_name == "point4"):
        skiprows = 258
    elif (ex_name == "point5"):
        skiprows = 338
    else:
        print("Not found the correspondent 'skiprows' int\n")
        sys.exit()
    
    traj = []
    
    for i in range(1, n_Conf + 1):
        coords_file = coords_path + "/Conf" + str(i).zfill(3) + ".xyz"
        traj.append(np.loadtxt(coords_file, skiprows = skiprows))
        
    new_file = coords_path + "/Adatom_traj.txt"
    np.savetxt(new_file, traj)
    
    return new_file