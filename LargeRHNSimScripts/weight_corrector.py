# -*- coding: utf-8 -*-
"""
This script corrects the weights obtained in simulation
"""
import sys
import numpy as np
import awkward as ak
import pickle 

sys.path.insert(0,'../../MATHUSLA_FastSim')
import Helpers.functions as hf

sys.path.insert(0, '../../FastSim_Additions')
from Additions import initiate_detector, get_ctaus

def uniform_error_correction(momentum, mass, ctau, detector_benchmark):
    multiplicative_correction = hf.get_weight(momentum, mass, ctau, detector_benchmark)[0]
    
    return multiplicative_correction

def main(filename):
    
    detector_benchmark = initiate_detector()
    
    with open(filename, 'rb') as f:
        data = pickle.load(f)['Data']
        
    flavour, _, mass, mixing_sq = filename[:-7].split('_')[-4:]
    mass = float(mass)
    mixing_sq = float(mixing_sq)
    ctau = get_ctaus(mass, mixing_sq, f'../../AllSimData/MATHUSLA_LLPfiles_RHN_{flavour}/RHNctau{flavour}.dat')
    
    four_p = data['momentum'] 
    
    weight_correction = np.array([uniform_error_correction(mom, mass, ctau, detector_benchmark)
                         for mom in four_p])
    
    data['weight'] = data['weight'] * weight_correction 

    return data['weight']
if __name__ == '__main__':
    filename = 'Finished_Sim_Ue/sim_Ue_Bmeson_1.1283_1.9306977288832498e-06.pickle'
    result = main(filename)
