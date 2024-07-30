#note: running this script might take a while, as we run both the vanzadelhoff 1a and 1b models 4 times each.
import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../tests/data/'
moddir = f'{curdir}/../tests/models/'
resdir = f'{curdir}/../tests/results/'

import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte
import matplotlib
from matplotlib.lines import Line2D

matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['axes.labelsize'] =  "large"
matplotlib.rcParams['legend.fontsize'] = "large"

from scipy.interpolate import interp1d


#Note: parameters correspond with the ones used for benchmarking the vanzadelhoff 1a/b models in 1D
dimension = 1
npoints   = 250
nrays     = 200
nspecs    = 5
nlspecs   = 1
nquads    = 11

r_in   = 1.0E13   # [m]
r_out  = 7.8E16   # [m]
nH2_in = 2.0E13   # [m^-3]
temp   =  20.00   # [K]
turb   = 150.00   # [.]
npops = 2

get_X_mol = {
    'a' : 1.0E-8,
    'b' : 1.0E-6
}

rs = np.logspace (np.log10(r_in), np.log10(r_out), npoints, endpoint=True)

def plot():

    modelName = f'vanZadelhoff_1a_1D'

    mem_lims = [4,8,16,32]
    fig, ax = plt.subplots(2,1)
    plt.sca(ax[0])

    for mem_lim in mem_lims:
        rel_change_max_without = np.load(f"ng_accel_comparison_no_adaptive_default_{modelName}_memlim_{mem_lim}.npy")
        print(rel_change_max_without)

        #without adaptive ng-accel
        plt.plot(1+np.arange(len(rel_change_max_without)), rel_change_max_without, linestyle="dashed")
        # plt.plot(1+np.arange(len(rel_change_max_without)), rel_change_mean_without, label="mean")
        #plt.xlabel("iteration")
        plt.ylabel("relative change")
        plt.legend()
        plt.yscale("log")
        
    #reset plot colors
    plt.gca().set_prop_cycle(None)

    for mem_lim in mem_lims:
        rel_change_max_with = np.load(f"ng_accel_comparison_adaptive_default_{modelName}_memlim_{mem_lim}.npy")
        print(rel_change_max_with)

        plt.plot(1+np.arange(len(rel_change_max_with)), rel_change_max_with, label=f'Nsteps, Nmax = ${mem_lim}$')
        # plt.plot(1+np.arange(len(rel_change_max_with)), rel_change_mean_with, label="mean, adaptive", linestyle="dashed")
        #plt.xlabel("Iteration")
        plt.ylabel("Max relative change")
        plt.legend(loc = "upper right")
        plt.yscale("log")
    
    plt.title("Problem 1a")
    
    modelName = f'vanZadelhoff_1b_1D'
    plt.sca(ax[1])
    plt.title("Problem 1b")
    
    for mem_lim in mem_lims:
        rel_change_max_without = np.load(f"ng_accel_comparison_no_adaptive_default_{modelName}_memlim_{mem_lim}.npy")
        print(rel_change_max_without)

        #without adaptive ng-accel
        plt.plot(1+np.arange(len(rel_change_max_without)), rel_change_max_without, linestyle="dashed")
        # plt.plot(1+np.arange(len(rel_change_max_without)), rel_change_mean_without, label="mean")
        plt.xlabel("iteration")
        plt.ylabel("relative change")
        plt.legend()
        plt.yscale("log")
        
    #reset plot colors
    plt.gca().set_prop_cycle(None)

    for mem_lim in mem_lims:
        rel_change_max_with = np.load(f"ng_accel_comparison_adaptive_default_{modelName}_memlim_{mem_lim}.npy")
        print(rel_change_max_with)

        plt.plot(1+np.arange(len(rel_change_max_with)), rel_change_max_with, label=f'Nsteps = ${mem_lim}$')
        # plt.plot(1+np.arange(len(rel_change_max_with)), rel_change_mean_with, label="mean, adaptive", linestyle="dashed")
        plt.xlabel("Iteration")
        plt.ylabel("Max relative change")
        
        handles = [Line2D([0], [0], label='Classical Ng', linestyle = "dashed", color="black"), Line2D([0], [0], label='Adaptive Ng', linestyle = "solid", color="black")]
        plt.legend(handles=handles, loc = "upper right")
        plt.yscale("log")
        
    plt.tight_layout()
    
    plt.savefig(f'ng_accel_comparison_adaptive_default_both_vanzadelhoff_both.pdf')
    plt.show()

    #plt.savefig(f'ng_accel_comparison_adaptive_default_{modelName}_multip le_memlims_remove_n_its_{remove_N_its}-{timestamp}.png')
    #plt.show() #will interrupt the program flow, so commented out

def run_test (nosave=True):

#I suggest uncommenting either the 'a' part or the 'b' part, as plt.show() will interrupt the python script
    plot()

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
