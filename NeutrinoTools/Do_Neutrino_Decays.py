# -*- coding: utf-8 -*-
"""
Intakes a neutrino file and outputs the tailored event data
"""
import sys
import numpy as np
import awkward as ak
from scipy.interpolate import interp1d
import uproot 
import pickle
import timeit

sys.path.insert(0, '../../MATHUSLA_FastSim/')
import DetectorSimulation.Detector as Detector
import DetectorSimulation.Particle as Particle

sys.path.insert(0,'../../FastSim_Additions')
from Additions import initiate_detector

import nu_data_structure as nds
from Repickling import main as repickle

#####Do Not Use, Doesn't Work #####
def get_event_weight(Ev, weight_file):
    '''
    Get weight of events based off of the incoming neutrino energy

    Parameters
    ----------
    Ev : np.ndarray
        Array of incoming energies 

    Returns
    -------
    np.ndarray
    '''
    weights = repickle(weight_file)['weights']
    energies = np.arange(0, len(weights), 1) / 4
    interp_func = interp1d(energies, weights, kind = 'previous')
    
    return interp_func(Ev)
    
def get_event_positions(number_of_events):
    '''
    Randomly generates event positions within MATHUSLA decay volume

    Parameters
    ----------
    number_of_events : int
        Number of positions to generate

    Returns
    -------
    np.ndarray
        An Nx3 array containing all vertex positions

    '''
    detector_benchmark = initiate_detector()
    
    position_fractions = np.random.rand(number_of_events, 3)
    
    scale_corrections = np.array([detector_benchmark.config.decay_x_dim,
                                  detector_benchmark.config.decay_y_dim,
                                  detector_benchmark.config.decay_z_max - detector_benchmark.config.decay_z_min])
    
    shift_corrections = np.array([detector_benchmark.config.decay_x_min, 
                                  detector_benchmark.config.decay_y_min, 
                                  detector_benchmark.config.decay_z_min])
    
    return position_fractions * scale_corrections + shift_corrections

def format_detector_particle(daughter):
    return [*daughter.particle.momentum,
            daughter.particle.pid,
            len(daughter.tracker_hits),
            daughter._visibility]

def do_decays(Ev, weights, pos, daughters):
    start = timeit.default_timer()
    
    data = nds.initiate_structure()
    detector_benchmark = initiate_detector()
    
    for i, event in enumerate(daughters):
        
        per_event = {'Ev':Ev[i], 'weight':weights[i], 'position':pos[i], 'daughters':[]}
        
        for daughter in event:
            detector_benchmark.clear_detector()
            
            particle = Particle.Particle(pos[i], daughter[:4], daughter[4])
            detector_benchmark.new_particle_event(particle)
            formatted = format_detector_particle(detector_benchmark.return_current_particle())
            if i == 0:
                print('Particle Successfully Decayed')
            per_event['daughters'].append(formatted)
            
        data = nds.append_new_event(per_event, data)
        
        if i % 10000 == 0:
            print(f'Event {i} Appended Successfully')
    end = timeit.default_timer()
    print(f'Decay time: {end - start}')
    return data
            
def main(genie_file, weight_file, outfile):
    #Load in entire data set
    start = timeit.default_timer()
    genie_events = nds.pull_from_root(genie_file)
    end = timeit.default_timer()
    print(f'Loaded ROOT Successfully (time: {end - start})')
    
    #Extract important elements
    #event_weights = get_event_weight(genie_events['Ev'], weight_file)
    event_weights = repickle(weight_file)['weights']
    print('Calculated Event Weights Successfully')
    event_positions = get_event_positions(len(genie_events['iev']))
    print('Generated Positions Successfully')
    daughters = nds.transpose_daughters(genie_events['Ef'], genie_events['pxf'], genie_events['pyf'],
                                        genie_events['pxf'], genie_events['pdgf'])
    
    all_data = do_decays(genie_events['Ev'], event_weights, event_positions, daughters)
    
    all_data = nds.freeze_structure(all_data)
    
    with open(outfile, 'wb') as f:
        
        pickle.dump(all_data,f)
        
    return all_data
    
if __name__ == '__main__':
    
    param_tuples = []
    for flavour in ['nue', 'nuebar', 'numu', 'numubar']:
        param_tuples.append((flavour,0,1))
        param_tuples.append((flavour,1,2))
        param_tuples.append((flavour,2,10))
        
    for params in param_tuples:
        genie_file = f'../../AllSimData/neutrino_files/{params[0]}_{params[1]}-{params[2]}Gev_100kev_gntp.0.gst.root'
        weight_file = f'../../AllSimData/neutrino_files/{params[0]}_{params[1]}-{params[2]}Gev_100kev_gntp.0.gst.pickle'
        outfile = f'../NeutrinoResults/{params[0]}_{params[1]}-{params[2]}_results.pickle'
        print(params) 
        main(genie_file, weight_file, outfile)
        
    print('Process Complete')
