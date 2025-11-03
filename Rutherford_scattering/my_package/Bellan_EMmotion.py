# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 14:16:38 2025

@author: SS882604
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['axes.labelsize'] = 13
mpl.rcParams['axes.titlesize'] = 13
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['legend.fontsize'] = 13
mpl.rcParams['figure.figsize'] = (8, 8)
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.bbox'] = 'tight'

#%% --------------- SIMULATOR STUFF ---------------    

# Core resolver

def compute_Vel (Vel, Omega, Sigma, dt):
    A = Omega * dt /2
    C = Vel + dt * (Sigma + np.cross(Vel, Omega) / 2)
    New_Vel = (C + A * np.dot(C, A) - np.cross(A, C) ) / (1 + np.dot(A, A))
    
    return New_Vel

# Used for the Rutherford scattering

def compute_E (Pos, q_fixed):
    r_2 = np.dot(Pos, Pos)
    r = np.sqrt(r_2)
    
    E_tot = q_fixed / (4 * cst.pi * cst.epsilon_0 * r_2)
    E_x = E_tot * Pos[0] / r
    E_y = E_tot * Pos[1] / r
    E_z = 0 
    
    return np.array([E_x, E_y, E_z])

def MRUA (time, acc):
    return 0.5 * acc * time**2

def HM (time, radius, w):
    return radius * np.sin(w * time)

def ExB (time, vel, w):
    # See: https://physicspages.com/pdf/Electrodynamics/Cycloid%20motion%20in%20crossed%20electric%20and%20magnetic%20fields.pdf
    return vel * time - vel/w * np.sin(w * time)
    
def nothing (time):
    return np.zeros_like(time)

#%% --------------- MATRIX OPERATIONS ---------------

def give_RotMat (axis, angle):
    
    match axis:
        case 'x':
            mat = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
            return mat
        case 'y':
            mat = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
            return mat
        case 'z':
            mat = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
            return mat
        case _:
            print("\n------------------------------------")
            print("      AXIS OF ROTATION NOT GIVEN     ")
            print("------------------------------------\n")
            
            mat = np.array([0, 0, 0], [0, 0, 0], [0, 0, 0])
            return mat
        
def compute_MatVec (mat, vec):
    a = np.zeros_like(vec)
    
    for i in range(0, len(a)):
        for j in range(0, len(a)):
            a[i] += mat[i,j]*vec[j]
            
    return a

#%% --------------- GRAPHICAL VISUALIZATIONS ---------------

# Obsolete function: works only if init_Vel = [..., 0, 0] and init_Pos = [..., ..., 0]

def plot_2DRutherford (all_Pos, all_cmPos1, all_cmPos2, all_Pos1, all_Pos2, cm_Vel):
    
    fig, [axE, axCM, axLAB] = plt.subplots(3, 1)

    axE.plot(all_Pos[:, 0], all_Pos[:, 1], '.b')
    axE.scatter(0, 0, marker="x", color="r")
    axE.axis('equal')
    axE.set_xlabel("x position [m]")
    axE.set_ylabel("y position [m]") 
    axE.set_title("Equivalent Trajectory")
    # axE.set_xlim(-50, +50)
    # axE.set_ylim(-50, +50)
    axE.grid()

    axCM.plot(all_cmPos1[:, 0], all_cmPos1[:, 1], color="blue", label='First')
    axCM.plot(all_cmPos2[:, 0], all_cmPos2[:, 1], color="red", label='Second')
    axCM.annotate("", xytext=all_cmPos1[-4, 0:2], xy=all_cmPos1[-3, 0:2], 
                  arrowprops=dict(arrowstyle="-|>", color='blue'))
    axCM.annotate("", xytext=all_cmPos2[-4, 0:2], xy=all_cmPos2[-3, 0:2], 
                  arrowprops=dict(arrowstyle="-|>", color='red'))
    axCM.set_xlabel("x position [m]")
    axCM.set_ylabel("y position [m]")
    axCM.set_title("CM Trajectories")
    axCM.grid()
    axCM.legend()

    axLAB.plot(all_Pos1[:, 0], all_Pos1[:, 1], color="blue", label='First - ')
    axLAB.plot(all_Pos2[:, 0], all_Pos2[:, 1], color="red", label='Second - ')
    axLAB.annotate("", xytext=all_Pos1[-4, 0:2], xy=all_Pos1[-3, 0:2], 
                  arrowprops=dict(arrowstyle="-|>", color='blue'))
    axLAB.annotate("", xytext=all_Pos2[-4, 0:2], xy=all_Pos2[-3, 0:2], 
                  arrowprops=dict(arrowstyle="-|>", color='red'))
    axLAB.set_xlabel("x position [m]")
    axLAB.set_ylabel("y position [m]")
    axLAB.set_title(r"LAB Trajectories: $v_{CM}$ = " + f"[{cm_Vel[0]:.2e}, {cm_Vel[1]:.2e}, {cm_Vel[2]:.2e}]")
    axLAB.grid()
    axLAB.legend()

    plt.tight_layout()
    plt.show()
    
# General plotting function
    
def plot_3DRutherford(all_Pos1, all_Pos2, title):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    ax.plot(all_Pos1[:,0], all_Pos1[:,1], all_Pos1[:,2], c='Blue', label='First')
    ax.plot(all_Pos2[:,0], all_Pos2[:,1], all_Pos2[:,2], c='Red', label='Second')
    
    ax.scatter(all_Pos1[-1,0], all_Pos1[-1,1], all_Pos1[-1,2], marker='o', color='Blue', s=64)
    ax.scatter(all_Pos2[-1,0], all_Pos2[-1,1], all_Pos2[-1,2], marker='o', color='Red', s=64)
    
    ax.legend()
    ax.set_title(title)
    ax.set_xlabel('X coord [m]')
    ax.set_ylabel('Y coord [m]')
    ax.set_zlabel('Z coord [m]')
    
    plt.show()
    
# Plot the simulation results step by step    
    
def plot_cmFrame (all_Pos):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    
    ax.scatter(all_Pos[:,0], all_Pos[:,1], all_Pos[:,2], c='Blue')
    ax.scatter(0, 0, 0, marker='x', c='Red')
    
    ax.set_title("EQUIVALENT problem - Trajectory")
    ax.set_xlabel('X coord [m]')
    ax.set_ylabel('Y coord [m]')
    ax.set_zlabel('Z coord [m]')
    
    plt.show()
    
# Plot a vector in 3D space, in a 'fancy' way
    
def plot_Vec(fig, ax, a, color, label):
    ax.plot([0, a[0]], [0, a[1]], [0, a[2]], color = color, label = label)
    ax.scatter(a[0], a[1], a[2], marker="o", color = color, s = 100)
    ax.scatter(0, 0 , 0, marker="o", color='k', s=100)