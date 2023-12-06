# -*- coding: utf-8 -*-
"""
This script will make the plot of mixing angle vs the Probability of decay within
MATHUSLA.  

Potential parameters to change:
    - |Ue^2| range
    - Number of vectors per mixing angle
    - Mass of the B-meson
    - Which angle is being mixed
"""
import sys
import os

import numpy as np

from MixingAngle_vs_Yield_Plot import get_ctaus, decay_in_MATHUSLA, make_plot
from multiprocessing import Pool

sys.path.insert(0, '../MATHUSLA_FastSim/')
import DetectorSimulation.Detector as Detector
import DetectorSimulation.llp_gun_new as lg
from Helpers.functions import *

sys.path.insert(0,'../FastSim_Additions/')
from Additions import initiate_detector

file = '../SimulationData/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv'
fv_path = os.path.join(os.getcwd(), file)
mass = float(file.split('_')[-1][:-4])
mixing_sq = np.logspace(-0.5,-7, 100) 
lengths = get_ctaus(mass = mass, mixing_sq = mixing_sq)

detector_benchmark = initiate_detector()
phi_min, phi_max, theta_min, theta_max = get_detector_angles(detector_benchmark)

def multi_pros(lengths):
    
    p = Pool(3)
    all_vals = np.array(p.map(single_mix_sim,lengths))
    
    return all_vals.T[0], all_vals.T[1]
    

def single_mix_sim(ctau):
    in_detector = 0
    decay_in_detector = 0
    vectors = read_vectors(fv_path, 6000)
    for vec in vectors:
        four_p = vec[1:]
        llp_theta = get_theta(four_p[1:])
    
        if (llp_theta < theta_max) and (llp_theta > theta_min):
            #print('Particle in MATHUSLA')
            in_detector += 1
            detector_benchmark.clear_detector()
        
            rot_four_p = deal_with_phi(four_p, phi_min, phi_max)
        
            pack = get_weight(rot_four_p, mass, ctau, detector_benchmark)
            
            if pack is not None:
                if decay_in_MATHUSLA(detector_benchmark, pack[1]):
                    decay_in_detector += pack[0]
    return decay_in_detector, in_detector 

if __name__ == '__main__':
    
    num_decays, num_events = multi_pros(lengths)
    
    make_plot(mixing_sq, num_decays, num_events)
    
    
    
    


            

