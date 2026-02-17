# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 09:25:42 2026

@author: SS882604
"""

full_path = "Y:/My Drive/MIA_UNI/Magistrale/CMS/Working_folder/KMC/Ag epitaxy/"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.constants as cst

import os
import sys

mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['figure.titlesize'] = 20
mpl.rcParams['figure.titleweight'] = 'bold'
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['figure.figsize'] = (8, 6)
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.bbox'] = 'tight'

#%% FUNCTIONS

def computeSnap_theta_rough (snap_name):
    
    # Each row of "snap_name" stores the position of one atom (x, y, z)
    # In our case we're interested in only the z axis, thus the elevation of each atom
    # Therefore we extract only the 2-th column
    
    heights = np.loadtxt(full_path + folder_name + "/" + snap_name, dtype = np.int64,
                         skiprows = 2, usecols = 2)
    
    # Split 'heights' into bins and remove the multiple countings, because:
    #   the 1st bin contains the '# of columns with an atom at z = 0'
    #   the 2nd bin contains the '# of columns with an atom at z = 1'
    #   and so on...
    #
    # For example:
    # if max_h = 2, then len(h_distrib) = 3, thus we correct till the 1-th entry
    # of 'h_distrib'.
    #
    # Clearly, the number of columns of height 0 is given by '1st bin' - '2nd bin'
    # Everything works fine if after cleaning: np.sum(h_distrib) = 3600, that is LX*LY
    
    h_distrib = np.bincount(heights)   # Bins start from z = 0 and reach z = max_h
    max_h = len(h_distrib) - 1
    h_values = np.linspace(0, max_h, max_h + 1)   # Height of the corresponding bin
    
    # Removing multiple countings
    for i in range(0, max_h):
        h_distrib[i] = h_distrib[i] - h_distrib[i+1]
        
    print(np.sum(h_distrib))
        
    theta = np.sum(h_values*h_distrib) / np.sum(h_distrib)   # Experimental coverage
    rough = np.sqrt( np.sum( h_distrib*(h_values - theta)**2 ) / np.sum(h_distrib) )   # Experimental roughness
    
    return theta, rough



def stat_test(data):
    
    mean = np.mean(data)
    std = np.std(data)
    
    print(f"<...> = {mean:.2e} " + r"+-" + f" {std:.2e}")
    
    return mean, std
    


def plot_data (fig, ax, plot_type, label, ref = True, ref_label = "Expected", phi = 0.2):
    
    data = np.loadtxt(full_path + folder_name + "/Output.txt")
    time, THETA, tau, theta, rough = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    
    match plot_type:
        case "tau_d":
            mean_tau = np.mean(tau)
            tau_hist, bins = np.histogram(tau, bins = 40, range = (0, 0.015))
            tau_axis = np.delete(bins, -1, axis = 0) + (bins[1] - bins[0])/2
            d_tau = tau_axis[1] - tau_axis[0]
            
            # ----- Generate the exponential distribution by extractions
            
            # if (ref == True):
            #     k_depo = 0.2*60*60
            #     scale = 1e1
            #     RNG = np.random.default_rng(seed = 12345678900987654321)
            #     temp = np.ones( int(np.size(tau) * scale) ) * 1/k_depo
            #     exponential = RNG.exponential(temp)
            #     exp_hist, bins = np.histogram(exponential, bins = bins)
            #     ax.stairs(exp_hist/scale, bins, label = "Expected", color = "k",
            #               fill = True, alpha = 0.05)
                
            # ----- Generate the exponential distribution using theory
            
            if (ref == True):
                k_depo = phi * 60 * 60
                x_axis = np.linspace(0, np.max(tau), 10000)
                exp_hist = np.zeros(np.size(bins) - 1)

                # Obtain the area covered by the distribution between bins[i+1] and bins[i]
                for i in range(0, np.size(exp_hist)):
                    exp_hist[i] = np.exp(- k_depo * bins[i]) - np.exp(- k_depo * bins[i+1])

                exp_hist *= np.size(tau)   # <--- From probability to number of samples
                
                ax.plot(tau_axis, exp_hist, 'r-', label = ref_label)
            
            ax.bar(tau_axis, tau_hist, width = d_tau, edgecolor = 'k', alpha = 0.5, label = label)
            # ax.stairs(tau_hist, bins, label = label, linewidth = 1.7)   # <--- Stair plot: BAD LOOKING
            # ax.plot(tau_axis, tau_hist, '-', label = label)   # <--- Scatter plot: BAD LOOKING
            
            ax.set_title("Escape time distribution")
            # + "\n" + r"$<\tau> = $" + f"{mean_tau:.2e} [s]  vs  " + 
            #             r"${k_{DEPO}}^{-1} = $" + f"{1/k_depo:.2e} [s]"
            ax.set_xlabel(r"$\tau$ values [s]")
            ax.set_ylabel("# of occurrences")
            ax.set_yscale('log')
            ax.ticklabel_format(axis = 'x', style = 'sci', scilimits = (0, 0))
            ax.set_facecolor("lightgray")
            ax.grid(True, color = 'gray')
            
            # ax.text(x = 0.008, y = 3000, s = r"$<\tau> = $" + f"{mean_tau:.2e}\n" 
            #                                   + r"$\frac{1}{k_{DEPO}} = $" + f"{1/(0.2*3600):.2e}",
            #         fontsize = 'large', bbox = dict(boxstyle='round', facecolor='white', edgecolor = 'k'))
            
            # ----- Perform rudimentary statistichal tests
            
            stat_test(tau)
            
        case "tau_d_variable_bins":
            mean_tau = np.mean(tau)
            tau_hist, bins = np.histogram(tau, bins = 40)
            tau_axis = np.delete(bins, -1, axis = 0) + (bins[1] - bins[0])/2
            d_tau = tau_axis[1] - tau_axis[0]
            
            # ----- Generate the exponential distribution by extractions
            
            # if (ref == True):
            #     k_depo = 0.2*60*60
            #     scale = 1e1
            #     RNG = np.random.default_rng(seed = 12345678900987654321)
            #     temp = np.ones( int(np.size(tau) * scale) ) * 1/k_depo
            #     exponential = RNG.exponential(temp)
            #     exp_hist, bins = np.histogram(exponential, bins = bins)
            #     ax.stairs(exp_hist/scale, bins, label = "Expected", color = "k",
            #               fill = True, alpha = 0.05)
                
            # ----- Generate the exponential distribution using theory
            
            if (ref == True):
                k_depo = phi * 60 * 60
                x_axis = np.linspace(0, np.max(tau), 10000)
                exp_hist = np.zeros(np.size(bins) - 1)

                # Obtain the area covered by the distribution between bins[i+1] and bins[i]
                for i in range(0, np.size(exp_hist)):
                    exp_hist[i] = np.exp(- k_depo * bins[i]) - np.exp(- k_depo * bins[i+1])

                exp_hist *= np.size(tau)   # <--- From probability to number of samples
                
                ax.plot(tau_axis, exp_hist, 'r-', label = ref_label)
            
            ax.bar(tau_axis, tau_hist, width = d_tau, edgecolor = 'k', alpha = 0.5, label = label)
            # ax.stairs(tau_hist, bins, label = label, linewidth = 1.7)   # <--- Stair plot: BAD LOOKING
            # ax.plot(tau_axis, tau_hist, '-', label = label)   # <--- Scatter plot: BAD LOOKING
            
            ax.set_title("Escape time distribution")
            # + "\n" + r"$<\tau> = $" + f"{mean_tau:.2e} [s]  vs  " + 
            #             r"${k_{DEPO}}^{-1} = $" + f"{1/k_depo:.2e} [s]"
            ax.set_xlabel(r"$\tau$ values [s]")
            ax.set_ylabel("# of occurrences")
            ax.set_yscale('log')
            ax.ticklabel_format(axis = 'x', style = 'sci', scilimits = (0, 0))
            ax.set_facecolor("lightgray")
            ax.grid(True, color = 'gray')
            
            # ax.text(x = 0.008, y = 3000, s = r"$<\tau> = $" + f"{mean_tau:.2e}\n" 
            #                                   + r"$\frac{1}{k_{DEPO}} = $" + f"{1/(0.2*3600):.2e}",
            #         fontsize = 'large', bbox = dict(boxstyle='round', facecolor='white', edgecolor = 'k'))
            
            # ----- Perform rudimentary statistichal tests
            
            stat_test(tau)
            
        case "tau_t_cov":
            ax.plot(theta, tau, '.--', label = label, linewidth = 0.5)
            ax.set_title(r"Exctracted values of $\tau$")
            ax.set_xlabel("Measured coverage [ML]")
            ax.set_ylabel("Escape time [s]")
            ax.grid(True)
            
        case "tau_t_tim":
            ax.plot(time, tau, '.--', label = label, linewidth = 0.5)
            ax.set_title(r"Exctracted values of $\tau$")
            ax.set_xlabel("Simulation time [s]")
            ax.set_ylabel("Escape time [s]")
            ax.grid(True)
            
        case "rough":
            
            if (ref == True):
                r = lambda theta: np.sqrt(theta * (1 - 1/3600)) 
                ax.plot(theta, r(theta), 'k', linewidth = 10, alpha = 0.2, label = ref_label)
            
            ax.plot(theta, rough, '-', label = label)
            ax.set_title("Roughness evolution")
            ax.set_xlabel("Measured coverage [ML]")
            ax.set_ylabel("Roughness [ML]")
            ax.grid(True)
            
        case "theta":
            
            if (ref == True):
                t = lambda THETA: THETA
                ax.plot(THETA, t(THETA), 'k',  linewidth = 10, alpha = 0.2, label = ref_label)
            
            ax.plot(THETA, theta, '-', label = label)
            ax.set_title("Coverage comparison")
            ax.set_xlabel("Expected coverage [ML]")
            ax.set_ylabel("Measured coverage [ML]")
            ax.grid(True)
            
    return fig, ax
    
def tau_oscillations (folder_name, label):
    data = np.loadtxt(full_path + folder_name + "/Output.txt")
    time, tau, theta, rough, k_depo, k_diff = data[:, 0], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True)
    
    fig.suptitle("Tau oscillations")
    
    ax1.plot(time, tau, 'g.', label = label, linewidth = 0.5)
    ax1.set_title(r"Exctracted values of $\tau$")
    ax1.set_ylabel(r"$\tau$ [s]")
    ax1.legend()
    ax1.grid(True)
    
    ax2.plot(time, k_diff, 'b.', label = r"Total $k_{DIFF}$")
    ax2.set_title("Transition rate time evolution")
    ax2.plot(time, k_depo, 'r.', label = r"$k_{DEPO}$")
    ax2.set_ylabel("Rate [1/s]")
    ax2.legend()
    ax2.grid(True)
    
    ax3.plot(time, rough, 'g-', label = label)
    ax3.set_title("Roughness time evolution")
    ax3.set_ylabel("Rough [ML]")
    ax3.set_xlabel("Simulation time [s]")
    ax3.legend()
    ax3.grid(True)
    
    fig.show()

#%% 1st PART ANALYSIS

# folder_name = "1st_part/02_Feb03_23_46_phi_0.2"

sims = ["1st_part/02_Feb05_16_47_phi_0.2",
        "1st_part/02_Feb05_16_48_phi_0.2",
        "1st_part/02_Feb05_16_49_phi_0.2",
        "1st_part/02_Feb05_16_50_phi_0.2",
        "1st_part/02_Feb05_16_51_phi_0.2"]

references = [True, 
              False,
              False,
              False,
              False]

labels = ["Seed #1",
          "Seed #2",
          "Seed #3",
          "Seed #4",
          "Seed #5"]

'''Analyze the tau distribution - Morfology of the deposition (cool Ovito visuals) - Rough vs Actual Cov || Expected Cov'''


# ----- ANALIZZA UN SINGOLO SNAP SALVATO

# theta, rough = computeSnap_theta_rough("SNAP0038.xyz")


# ----- GRAFICI

# fig, ax = plt.subplots()

# for (sim, label, ref) in zip(sims, labels, references):
    
#     folder_name = sim

#     fig, ax = plot_data(fig, ax, "tau_t", label = label, ref = ref, ref_label = "Expected")

# ax.legend()
# fig.show()


# ----- FUFFA

# fig, ax = plt.subplots()

# for (sim, label, i) in zip(sims, labels, [1, 2, 3, 4, 5]):
#     data = np.loadtxt(full_path + sim + "/Output.txt")
#     THETA, tau, theta, rough = data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    
#     mean, std = stat_test(tau)
    
#     ax.scatter(int(i), mean, label = label)
    
# ax.legend()
# ax.set_ylim(0.0012, 0.0015)
# ax.set_ylabel("Mean escape time [s]")
# ax.set_title("")
# ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0, 0))
# ax.grid()
# fig.show()
    

#%% 2nd PART ANALYSIS
sims = ["2nd_part/02_Feb06_16_33_phi_0.0001",
        "2nd_part/02_Feb06_16_26_phi_0.001",
        "2nd_part/02_Feb06_16_05_phi_0.01",
        "2nd_part/02_Feb06_15_11_phi_0.1",
        "2nd_part/02_Feb06_15_53_phi_1.0",
        "2nd_part/02_Feb06_15_59_phi_10.0"]

references = [False, 
              False,
              False,
              False,
              False,
              False]

labels = [r"$\phi$ = 0.0001 [ML/s]",
          r"$\phi$ = 0.001 [ML/s]", 
          r"$\phi$ = 0.01 [ML/s]",
          r"$\phi$ = 0.1 [ML/s]",
          r"$\phi$ = 1.0 [ML/s]",
          r"$\phi$ = 10 [ML/s]"]

phis = [0.0001, 0.001, 0.01, 0.1, 1.0, 10]

# ----- GRAFICI

fig, ax = plt.subplots()

for (sim, label, ref) in zip(sims, labels, references):
    
    folder_name = sim

    fig, ax = plot_data(fig, ax, "tau_d_variable_bins", label = label, ref = ref, ref_label = "Deposition only")

ax.legend()
fig.show()

# ----- SINGOLO GRAFICO

# fig, ax = plt.subplots()
# folder_name = sims[5]
# label = labels[5]
# ref = True
# phi = phis[5]
# fig, ax = plot_data(fig, ax, "tau_d_variable_bins", label = label, ref = ref, ref_label = "Deposition only", phi = phi)
# ax.legend()
# fig.show()

# ----- TAU OSCILLATIONS

# tau_oscillations(sims[0], labels[0])

#%% IN SEARCH OF NEW T

# ----- FORMULAS

# K = 1/11603
# k_depo = 1 * 3600

# def k_diff (n, T):
#     return 4 * 1e13 * np.exp( - (4+n) * 0.345 / (K*T) )

# k_diff(0, 764.3) / k_diff(4, 764.3)

#%% IN SEARCH OF NEW T

# ----- PLOTS: comparison between seeds

# sims = ["2nd_part/02_Feb06_16_05_phi_0.01",
#         "2nd_part/02_Feb08_11_23_phi_0.01",
#         "2nd_part/02_Feb08_11_31_phi_0.01"]

# sims = ["2nd_part/02_Feb06_18_17_phi_1.0",
#         "2nd_part/02_Feb08_11_25_phi_1.0",
#         "2nd_part/02_Feb08_11_33_phi_1.0"]

sims = ["2nd_part/02_Feb08_11_17_phi_1.0",
        "2nd_part/02_Feb08_11_28_phi_1.0",
        "2nd_part/02_Feb08_11_37_phi_1.0"
        ]

references = [False,
              False,
              False]

labels = ["Seed #1", "Seed #2", "Seed #3"]

fig, ax = plt.subplots()

for (sim, label, ref) in zip(sims, labels, references):
    
    folder_name = sim

    fig, ax = plot_data(fig, ax, "rough", label = label, ref = ref, ref_label = "Deposition only")

# ax.set_title("Roughness evolution\n" + r"$\phi$ = 0.01 [ML/s] - T = 650 [K]")
ax.set_title("Roughness evolution\n" + r"$\phi$ = 1.0 [ML/s] - T = 764 [K]")
ax.legend()
fig.show()

#%% IN SEARCH OF NEW T

# ----- PLOTS: comparison between temperatures within same seed

# sims = ["2nd_part/02_Feb06_16_05_phi_0.01",   # Seed #1
#         "2nd_part/02_Feb06_18_17_phi_1.0",
#         "2nd_part/02_Feb08_11_17_phi_1.0"]

# sims = ["2nd_part/02_Feb08_11_23_phi_0.01",   # Seed #2
#         "2nd_part/02_Feb08_11_25_phi_1.0",
#         "2nd_part/02_Feb08_11_28_phi_1.0"]

sims = ["2nd_part/02_Feb08_11_31_phi_0.01",   # Seed #3
        "2nd_part/02_Feb08_11_33_phi_1.0",
        "2nd_part/02_Feb08_11_37_phi_1.0"]

references = [False,
              False,
              False]

labels = [r"$\phi$ = 0.01 [ML/s] - T = 650 [K]",
          r"$\phi$ = 1.0 [ML/s] - T = 800 [K]",
          r"$\phi$ = 1.0 [ML/s] - T = 764 [K]"]

fig, ax = plt.subplots()

for (sim, label, ref) in zip(sims, labels, references):
    
    folder_name = sim

    fig, ax = plot_data(fig, ax, "rough", label = label, ref = ref, ref_label = "Deposition only")

ax.set_title("Roughness evolution: seed #3")
ax.legend()
fig.show()

#%% IN SEARCH OF NEW T

# ----- PLOTS: comparison between temperatures with same PHI

sims = ["2nd_part/02_Feb06_15_53_phi_1.0",   # Seed #1
        "2nd_part/02_Feb06_18_17_phi_1.0",
        "2nd_part/02_Feb08_11_17_phi_1.0"]

references = [False,
              False,
              False]

labels = ["T = 650 [K]",
          "T = 800 [K]",
          "T = 764 [K]"]

fig, ax = plt.subplots()

for (sim, label, ref) in zip(sims, labels, references):
    
    folder_name = sim

    fig, ax = plot_data(fig, ax, "rough", label = label, ref = ref, ref_label = "Deposition only")

ax.set_title(r"Roughness evolution: $\phi$ = 1.0 [ML/s]")
ax.legend()
fig.show()
