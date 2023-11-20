# -*- coding: utf-8 -*-
"""
This script is used to calculate the mixing angle of RHN model
"""
import sys
import os

import numpy as np
import awkward as ak
import matplotlib.pyplot as plt

import pickle

from Closest_Approach_Histogram import get_reconstructable

sys.path.insert(0,'../../FastSim_Additions/')
import Additions as ad
import Event_Data_Structure as eds

plt.style.use('science')
plt.rcParams.update({'legend.frameon':True})

def organize_sim_info(fv_file, length_file, product_file, param_card_file = '../../MATHUSLA_FastSim/param_card_CDR.txt'):
    '''
    Takes in a sim pack and outputs the list of required information to run a simulation

    Parameters
    ----------
    fv_file : str
        The file containing four-vectors of LLPs coming from CMS
    length_file : str
        A mass to proper lifetime conversion file
    product_file : str
        File containing possible decay products, in the Hadron Decay format
    param_card_file : str, optional
        Param card for MATHUSLA geometry. The default is '../../MATHUSLA_FastSim/param_card_CDR.txt'.
   num_particles: int, optional
        Number of particles taken from fv_file. The default is 20000. -1 uses all.  

    Returns
    -------
    list
        A list characterizing simulation parameters

    '''
    #Relative Paths to important Information
    fv_path = os.path.join(os.getcwd(), fv_file)
    len_path = os.path.join(os.getcwd(), length_file)
    prod_path = os.path.join(os.getcwd(), product_file)

    mass = float(fv_file.split('_')[-1][:-4])
    mixing = 0.004 #This will eventually need to be changed 
    ctau = ad.get_ctaus(mass = mass, mixing = mixing, length_file = len_path)
    detector_benchmark = ad.initiate_detector('../../MATHUSLA_FastSim/param_card_CDR.txt')
    
    #The last integer is the simulation size
    return [fv_path, prod_path, ctau, mass, detector_benchmark, 50000]


def isolate_usable_3p(data, hit_threshold):
    '''
    Gives events which have a defined opening angle

    Parameters
    ----------
    data : dict
        Data dictionary from Event_Data_Structure
    hit_threshold : int
        Number of hits required for a reconstructable particle

    Returns
    -------
    ak.Array
        Events with defined opening angles
    bool_slicing_none_events : ak.Array
        bool awkard array to slice out events without defined opening angle

    '''
    reconstructable_daughters = get_reconstructable(data, hit_threshold)
    
    #Mask makes events which are false None types
    array_with_none_events = ak.mask(data['daughters'][:,:,1:4], reconstructable_daughters)
    #This makes a bool array for which events have at least two reconstructable tracks, so have defined opening angle
    bool_slicing_none_events = ak.sum(~ak.is_none(array_with_none_events, axis = 1), axis = 1) >= 2
    
    return ak.drop_none(array_with_none_events[bool_slicing_none_events]), bool_slicing_none_events
    
def eliminate_records(array_with_records):
    '''Records are what is returned by some awkard built-ins, they are hard to 
    work with and hard to get rid of'''
    build = ak.ArrayBuilder()

    to_fill = []
    
    for layer in array_with_records:
        to_fill = []
    
        for entry in layer:
            to_fill.append(list(entry.to_list()))
            
        build.append(to_fill)
    
    return build.snapshot()
    
def calculate_opening_angle(data):
    '''
    Calculates the opening angle for all simulated events in data

    Parameters
    ----------
    data : dict
        Data dictionary from Event_Data_Structure

    Returns
    -------
    np.ndarray
        The maximum opening angle for each simulated event
    weight_bool : ak.Array
        An array that tracked which events had a defined opening angle.  Used to slice weights

    '''

    usable_3p, weight_bool = isolate_usable_3p(data)
    
    #Finds all possible opening angle calculation pairs
    #Maps n x m x 3 array to n x mC2 x 2 x 3 array
    all_combs = eliminate_records(ak.combinations(usable_3p, 2, axis = 1))
    
    #Maps n x mC2 x 2 x 3 array to n x mC2 array
    dot_of_combs = ak.sum(ak.prod(all_combs, axis = 2), axis = 2)
    #Mapes n x mC2 x 2 x 3 array to n x mC2 x 2 array
    norm_of_combs = np.sqrt(ak.sum(all_combs**2, axis = 3))
    
    #Maps n x mC2 and n x mC2 x 2 array to n x mC2 array
    all_theta = np.arccos(dot_of_combs / ak.prod(norm_of_combs, axis = 2))
    
    #Gives maximum opening angle
    return ak.max(all_theta, axis = 1).to_numpy(), weight_bool
    

def make_plot(angles, weights):
    labels = ['RHN', 'HXX']
    
    fig, ax = plt.subplots(1,1, figsize = (6,4))
    for angle, label, weight in zip(angles, labels, weights):
        print(len(angle))

        ax.hist(angle, bins = np.arange(0, np.pi, 0.1), histtype = 'step', label = label, weights = weight)
    
    ax.grid()
    ax.legend(title = 'Model')
    ax.set_xlabel('Opening Angle of Vertex')
    ax.set_ylabel('Counts')
    ax.set_title('Comparing Opening Angles of Different Models')
    plt.savefig('Opening_Angle.pdf')
    

def main(sim_packs):
    '''
    Simulates events and tracks the opening angle of each event.  Produces a histrogram of the data.  

    Parameters
    ----------
    sim_packs : list
        Simulation parameters

    Returns
    -------
    None.

    '''
    
    angles = []
    weights = []
    
    for i,pack in enumerate(sim_packs):
        
        sim_info = organize_sim_info(*pack)
        vertices =  ad.time_process(ad.do_sim, sim_info)
        
        data = eds.to_data_structure(vertices)
        eds.data_to_pickle(data, f'Theta_Simulation_Data_Pack{i}.pickle')
        
        single_model_theta, weight_bool = calculate_opening_angle(data)
        angles.append(single_model_theta)
        weights.append(data['weight'][weight_bool])
        
    make_plot(angles, weights)
    
    with open('Calculated_Opening_Angles.pickle','wb') as f:
        pickle.dump((angles, weights), f)
    
    return angles, weights
    
if __name__ == '__main__':
    data = main([['../../SimulationData/ForMixingPlot/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv',
          '../../SimulationData/RHNctauUe.dat', '../../SimulationData/vN_to_all_0.314228_hadron_style.txt'],])

