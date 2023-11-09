# -*- coding: utf-8 -*-
"""
Creates a Plot confirming conservation of 4 momenta (and gives insight to fp rounding)

Possible parameters to vary are the Particle Mass, Flavour Coupling, and Mixing Angle
Also possible to vary simulation size, but probably not necessary

"""
import sys
import os
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import timeit

plt.style.use('science')

sys.path.insert(0, '../../MATHUSLA_FastSim/')
import DetectorSimulation.Detector as Detector
import DetectorSimulation.llp_gun_new as lg
from Helpers.functions import *

sys.path.insert(0,'../../FastSim_Additions/')
from Additions import initiate_detector
from run_simulation import do_sim
import Event_Data_Structure as eds

sys.path.insert(0, '../MixingAnglePlots/')
from MixingAngle_vs_Yield_Plot import get_ctaus


def time_process(func, args):
    '''A function that times the input funciton
    It returns the output of the function that it timed
    '''
    start = timeit.default_timer()
    output = func(*args)
    end = timeit.default_timer()
    
    print(f'{str(func).split()[1]}: Runtime {end - start} s')
    
    return output

def reformat(vertices, pickle_name = 'Conservation_Data.pickle'):
    '''Reforms list of vertices to data structure and outputs to a pickle'''
    data = eds.to_data_structure(vertices)
    eds.data_to_pickle(data, pickle_name)
    
    return data

def momentum_differential(data):
    '''Computes the difference in 4-momentum before and after decay at a vertex'''
    momentum_before = data['momentum']
    momentum_after = ak.sum(data['daughters'][:,:,:4], axis = 1).to_numpy()
    
    return momentum_before - momentum_after
    
def make_plots(four_p_diff):
    '''Makes the plot of the four momentum difference'''
    key = [r'$E$', r'$p_x$', r'$p_y$', r'$p_z$']
    
    #Loops through E, px, py, pz
    for i,difference in enumerate(four_p_diff.T):
        bins = np.arange(-0.5,0.50001, 0.01) #The 0.0001 is to ensure 5 is included
        plt.hist(difference, bins = bins, histtype = 'step', label = key[i])
    
    plt.title('Conservation of 4-Momentum')
    plt.xlabel(r'$p^{\mu}_{intial} - p^{\mu}_{final}$ (GeV)')
    plt.ylabel('Counts')
    plt.legend(frameon = True)
    plt.grid()
    
    #To change the output name
    plt.savefig('MomentumConservationConfirmation.pdf')
    
    
def main(fv_file = '../../SimulationData/ForMixingPlot/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv',
         length_file = '../../SimulationData/RHNctauUe.dat',
         product_file = '../../SimulationData/vN_to_all_0.314228_new.txt',
         mixing = 0.004, particle_number = 20000):
    '''Main takes all parameters to run the simulation excluding the number of points'''
    
    fv_path = os.path.join(os.getcwd(), fv_file)
    len_path = os.path.join(os.getcwd(), length_file)
    prod_path = os.path.join(os.getcwd(), product_file)
    
    mass = float(fv_file.split('_')[-1][:-4])
    ctau = get_ctaus(mass = mass, mixing = mixing,length_file = len_path)
    
    detector_benchmark = initiate_detector('../../MATHUSLA_FastSim/param_card_CDR.txt')

    vertices = time_process(do_sim, [fv_path, prod_path, ctau, mass, detector_benchmark, particle_number])
    
    data = time_process(reformat,[vertices])
    
    four_p_diff = momentum_differential(data)
    
    make_plots(four_p_diff)
    
    print('Complete')
        
if __name__ == '__main__':
    main()
    
    
