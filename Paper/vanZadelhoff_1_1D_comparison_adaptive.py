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


def create_model (a_or_b):
    """
    Create a model file for the van Zadelhoff 1 benchmark in 1D.
    """

    modelName = f'vanZadelhoff_1{a_or_b}_1D'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'

    X_mol = get_X_mol[a_or_b]

    def nH2 (r):
        return nH2_in * np.power(r_in/r, 2.0)

    def nTT (r):
        return X_mol  * nH2(r)

    model = magritte.Model ()
    model.parameters.set_spherical_symmetry(True)
    model.parameters.set_model_name        (modelFile)
    model.parameters.set_dimension         (dimension)
    model.parameters.set_npoints           (npoints)
    model.parameters.set_nrays             (nrays)
    model.parameters.set_nspecs            (nspecs)
    model.parameters.set_nlspecs           (nlspecs)
    model.parameters.set_nquads            (nquads)

    model.geometry.points.position.set([[r, 0, 0] for r in rs])
    model.geometry.points.velocity.set(np.zeros((npoints, 3)))

    model.chemistry.species.abundance = [[     0.0, nTT(r), nH2(r),  0.0,      1.0] for r in rs]
    model.chemistry.species.symbol    =  ['dummy0', 'test',   'H2', 'e-', 'dummy1']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)
    model = setup.set_Delaunay_boundary       (model)
    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_rays_spherical_symmetry (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return


def run_model (a_or_b, nosave=False):

    modelName = f'vanZadelhoff_1{a_or_b}_1D'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()
    
    remove_N_its = 0#number of iterations to ignore before starting lambda iteration
    #does not seem to make a significant impact on the performance, so keep at 0
    max_its = 100

    mem_lims = [4,8,16,32]
    plt.figure()

    for mem_lim in mem_lims:
    
        #reset computation for without adaptive ng
        model.compute_LTE_level_populations   ()
        model.parameters.use_adaptive_Ng_acceleration = False
        model.parameters.pop_prec = 1.0e-14
        #set ng acceleration params
        model.parameters.Ng_acceleration_mem_limit = mem_lim
        model.parameters.Ng_acceleration_remove_N_its = remove_N_its
        model.parameters.adaptive_Ng_acceleration_min_order = 4

        timer3 = tools.Timer('running model')
        timer3.start()
        model.compute_level_populations (True, max_its)
        timer3.stop()

        pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), npops))
        abun = np.array(model.chemistry.species.abundance)[:,1]
        rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)
        rel_change_max_without = np.array(model.error_max)
        rel_change_mean_without = np.array(model.error_mean)
        print(rel_change_max_without)
        print(rel_change_mean_without)

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
        
        #with adaptive ng-accel
        #reset computation for with adaptive ng
        model.compute_LTE_level_populations   ()
        model.parameters.use_adaptive_Ng_acceleration = True
        #set ng acceleration params
        model.parameters.Ng_acceleration_mem_limit = mem_lim
        model.parameters.Ng_acceleration_remove_N_its = remove_N_its

        timer4 = tools.Timer('running model')
        timer4.start()
        model.compute_level_populations (True, max_its)
        timer4.stop()


        pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), npops))
        abun = np.array(model.chemistry.species.abundance)[:,1]
        rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)
        rel_change_max_with = np.array(model.error_max)
        rel_change_mean_with = np.array(model.error_mean)
        print(rel_change_max_with)
        print(rel_change_mean_with)

        plt.plot(1+np.arange(len(rel_change_max_with)), rel_change_max_with, label=f'Nsteps = ${mem_lim}$')
        # plt.plot(1+np.arange(len(rel_change_max_with)), rel_change_mean_with, label="mean, adaptive", linestyle="dashed")
        plt.xlabel("Iteration")
        plt.ylabel("Max relative change")
        plt.legend(loc = "lower left")
        plt.yscale("log")

    plt.savefig(f'ng_accel_comparison_adaptive_default_{modelName}_multiple_memlims_remove_n_its_{remove_N_its}-{timestamp}.png')
    #plt.show() #will interrupt the program flow, so commented out

def run_test (nosave=True):

#I suggest uncommenting either the 'a' part or the 'b' part, as plt.show() will interrupt the python script
    create_model ('a')
    run_model    ('a', nosave)

    create_model ('b')
    run_model    ('b', nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
