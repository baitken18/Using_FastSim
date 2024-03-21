# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:28:46 2023

@author: baitk
"""

import sys 
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

sys.path.insert(0, '../../MATHUSLA_FastSim/')
import DetectorSimulation.Detector as Detector
import DetectorSimulation.llp_gun_new as lg
from Helpers.functions import *

sys.path.insert(0,'../../FastSim_Additions/')
from Additions import initiate_detector

def get_ctaus(mass, mixing_sq, length_file = '../../SimulationData/RHNctauUe.dat'):
    '''
    Converts mass to ctau for a given mixing angle (input mixing angle squared)

    Parameters
    ----------
    mass : float
        Mass of the RHN in GeV
    mixing_sq : np.ndarray
        Squared Mixing Angles to use
    length_file : str, optional
        Conversion from mass to length for U = 1. The default is '../SimulationData/RHNctauUe.dat'.

    Returns
    -------
    lengths : TYPE
        DESCRIPTION.

    '''
    #Loads conversion file
    mass_to_length = np.loadtxt(length_file, skiprows = 1)
    
    #Sets up an interpolation function converting mass to length
    interp_mass = interp1d(mass_to_length.T[0], mass_to_length.T[1])
    
    #Computes for U=1 then corrects for chosen parameter
    length_U1 = float(interp_mass(mass))
    lengths = length_U1 / mixing_sq
    return lengths

def decay_in_MATHUSLA(detector_benchmark, p_pos):
    '''
    Checks if the particle decays in MATHUSLA 

    Parameters
    ----------
    detector_benchmark : Detector.Detector
        The detector being used
    p_pos : Position of the particle
        Position of the decay

    Returns
    -------
    bool
        If particle was in detector

    '''
    p_pos = np.array(p_pos)
    
    #Establishes boundaries of the detector
    lower_bnds = np.array([detector_benchmark.config.decay_x_min,
                           detector_benchmark.config.decay_y_min,
                           detector_benchmark.config.decay_z_min])
    upper_bnds = np.array([detector_benchmark.config.decay_x_max,
                           detector_benchmark.config.decay_y_max,
                           detector_benchmark.config.decay_z_max])
    
    #This statement is true iff particle within all bounds
    return np.sum(p_pos > lower_bnds) == 3 and np.sum(p_pos < upper_bnds) == 3

def do_sim(fv_path, lengths, mass, detector_benchmark, particle_number = 1000):
    '''
    Completes the entire simulation

    Parameters
    ----------
    fv_path : str
        Where to get particles from
    lengths : np.ndarray
        The ctau parameters
    mass : float
        Mass of the RHN
    detector_benchmark : Detector.Detector
        Detector being used

    Returns
    -------
    decay_in_detector : np.ndarray
        Number of particles that decay in the detector for each mixing angle
    in_detector : np.ndarray
        Number of particles that enter the detector (but don't neccesarily decay) for 
                                                     each mixing angle

    '''
    phi_min, phi_max, theta_min, theta_max = get_detector_angles(detector_benchmark)
    
    in_detector = np.zeros(len(lengths))
    decay_in_detector = np.zeros(len(lengths))
    
    #Loops through all lengths
    for i,ctau in enumerate(lengths):
        print('New Decay Length = ', ctau)
        
        #Reads in n different particles from the fv_path file
        vectors = read_vectors(fv_path, particle_number)
        
        #Loops through vectors and gets important quantities
        for vec in vectors:
            four_p = vec[1:]
            llp_theta = get_theta(four_p[1:])
            
            #True if the theta angle is within bounds of the detector
            if (llp_theta < theta_max) and (llp_theta > theta_min):
                in_detector[i] += 1
                
                detector_benchmark.clear_detector()
                
                #Rotates phi angle to be within the detector
                rot_four_p = deal_with_phi(four_p, phi_min, phi_max)
                
                #get_weight determines the decay position from a distribution
                #This is the pack[1] entry
                pack = get_weight(rot_four_p, mass, ctau, detector_benchmark)
                
                #Checks that the particle hits the detector (it should, but you never know)
                if pack is not None:
                    #Checks that the decay position is within MATHUSLA
                    if decay_in_MATHUSLA(detector_benchmark, pack[1]):
                        decay_in_detector[i] += pack[0]
                    '''
                    detector_benchmark.clear_detector()
                    
                    llp = lg.get_llp('leptonic2body', mass, pack[1], pack[2],[11,11])
                    
                    detector_benchmark.new_vertex_event(llp)
                    if bool(detector_benchmark.vertex_is_in_decay_volume()):
                        decay_in_detector[i] += 1
                    '''
                    
    return decay_in_detector, in_detector
    
def make_plot(mixing_sq, num_decays, num_events, filename = 'MixingAngle_decay_probability.png'):
    '''
    Makes a "histogram" of event rate vs mixing angle

    Parameters
    ----------
    mixing : np.ndarray
        The mixing angles used 
    num_decays : np.ndarray
        Number of decays
    num_events : np.ndarray
        Number of events

    Returns
    -------
    None.

    '''
    fig, ax = plt.subplots(1,1, figsize = (6,4), constrained_layout = True)
    ax.step(mixing_sq, num_decays / num_events)
    ax.set_title('Decay Probability in MATHUSLA')
    ax.set_xlabel(r'$|U_e|^2$')
    ax.set_ylabel(r'$\frac{N_{decay}}{N_{enter}}$')
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.grid()
    plt.savefig(filename)
    plt.show()
    plt.close()

def main(fv_file = '../../SimulationData/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv',
         length_file = '../../SimulationData/RHNctauUe.dat',
         image_file = 'MixingAngle_decay_probability.png'):
    #Particle Information
    fv_path = os.path.join(os.getcwd(), fv_file)
    mass = float(fv_file.split('_')[-1][:-4])
    mixing_sq = np.logspace(-1,-12, 30) 
    lengths = get_ctaus(mass = mass, mixing_sq = mixing_sq)
    
    #Detector Information and Boundaries
    detector_benchmark = initiate_detector()

    num_decays, num_events = do_sim(fv_path, lengths, mass, detector_benchmark)
    
    print(num_decays)
    print(num_events)

    make_plot(mixing_sq, num_decays, num_events, image_file)
    
if __name__ == '__main__':
    main()