# -*- coding: utf-8 -*-
"""
This is the script that makes the mixing angle plots
"""
import sys 
import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pickle

plt.style.use('science')

from MixingAngle_vs_Yield_Plot import get_ctaus, decay_in_MATHUSLA, do_sim

sys.path.insert(0, '../MATHUSLA_FastSim/')
import DetectorSimulation.Detector as Detector
import DetectorSimulation.llp_gun_new as lg
from Helpers.functions import *

sys.path.insert(0,'../FastSim_Additions/')
from Additions import initiate_detector

def calculate_decays(fv_file = '../SimulationData/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv',
         length_file = '../SimulationData/RHNctauUe.dat', particle_number = 1000):
    #Particle Information
    fv_path = os.path.join(os.getcwd(), fv_file)
    mass = float(fv_file.split('_')[-1][:-4])
    mixing_sq = np.logspace(-1,-12, 50) 
    lengths = get_ctaus(mass = mass, mixing_sq = mixing_sq)
    
    #Detector Information and Boundaries
    detector_benchmark = initiate_detector()

    num_decays, num_events = do_sim(fv_path, lengths, mass, detector_benchmark, particle_number)
    
    return mass, mixing_sq, num_decays, num_events
    
def format_plot(ax):
    ax.set_title('Decay Probability in MATHUSLA')
    ax.set_xlabel(r'$|U_e|^2$')
    ax.set_ylabel(r'$\frac{N_{decay}}{N_{enter}}$')
    ax.set_xscale('log')
    ax.grid()
    
def main(fv_files = ['../SimulationData/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv'],
         length_file = '../SimulationData/RHNctauUe.dat',
         image_file = 'Final_Mixing_Angle_to_Decay.pdf',
         particle_number = 800):
    
    fig, ax = plt.subplots(1,1, figsize = (6,4))
    format_plot(ax)
    
    for fv_file in fv_files:
        mass, mixing_sq, num_decays, num_events = calculate_decays(fv_file, length_file, particle_number)
        ax.step(mixing_sq, num_decays / num_events, label = f'{mass:.3g}')
        
    ax.legend(frameon = True, title = 'RHN Mass (GeV)')
    plt.savefig(image_file)
    
    
if __name__ == '__main__':
    #Replace fv_files with a list of whichever files you want
    fv_files = glob.glob('../SimulationData/ForMixingPlot/RHN_Ue*.csv')
    main(fv_files)