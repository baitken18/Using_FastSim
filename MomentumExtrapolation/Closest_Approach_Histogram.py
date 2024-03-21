import sys
import os

import numpy as np
import awkward as ak
import matplotlib.pyplot as plt

sys.path.insert(0,'../../FastSim_Additions/')
import Additions as ad
import Event_Data_Structure as eds

plt.style.use('science')
plt.rcParams.update({'legend.frameon':True})

def organize_sim_info(fv_file, length_file, product_file, 
                      param_card_file = '../../MATHUSLA_FastSim/param_card_CDR.txt', num_particles = 20000):
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
    mixing = 0.004
    ctau = ad.get_ctaus(mass = mass, mixing = mixing, length_file = len_path)
    detector_benchmark = ad.initiate_detector('../../MATHUSLA_FastSim/param_card_CDR.txt')
    
    #The last integer is the simulation size
    return [fv_path, prod_path, ctau, mass, detector_benchmark, num_particles]

def get_reconstructable(data, hit_threshold = 4):
    '''
    Checks which events are reconstructable

    Parameters
    ----------
    data : dict
        Dictionary containing event data
    hit_threshold : int, optional
        Number of layers a visible particle needs to hit. The default is 4.

    Returns
    -------
    ak.Array
        A boolean array that slices out daughters which aren't visible

    '''
    enough_hits = data['daughters'][:,:,5] >= hit_threshold
    not_invisible = data['daughters'][:,:,6] > -1
    
    return enough_hits & not_invisible

def separate_dot(x,y):
    '''Takes dot product of two matrices along axis = 1'''
    return np.sum(x * y, axis = 1)

def proj(x,y):
    '''Projects x onto y'''
    return separate_dot(x,y) * y.T / np.linalg.norm(y, axis = 1)**2

def get_nearest_approach(position, momentum, unitary = True):
    '''Calculates closest approach distance for a particle using perp distance formula'''
    if unitary:
        momentum = momentum / ak.sum(momentum**2)**0.5
        
    projection = position - proj(position, momentum).T
    return np.linalg.norm(projection, axis = 1)
    
def calculate_shortest_distance(data, track_thresh = 2):
    '''
    Calculates the shortest distance between a reconstructed track and the IP

    Parameters
    ----------
    data : dict
        Data structure used for all simulations

    Returns
    -------
    shortest_distance : np.ndarray
        Closest distance to IP for each run
    reconstructable_position : np.ndarray
        A boolean array which slices only events with visible decays

    '''
    #These are boolean masks which tell when events are reconstructable
    #Second line checks that the vertex has at least one reconstructed daughter
    reconstructable_daughters = get_reconstructable(data)
    reconstructable_position = ak.sum(reconstructable_daughters, axis = 1) >= track_thresh
    
    reconstructed_3p = ak.fill_none(ak.mask(data['daughters'][:,:,1:4], reconstructable_daughters),
                                    ak.Array([0,0,0]),axis = 1)
    
    if 'Ev' in data.keys():
        total_reconstructed_momentum_usable = ak.sum(reconstructed_3p[reconstructable_position], axis = 1).to_numpy()
        data['position'] = data['position'].reshape(len(data['weight']), 3)
    else:
        total_reconstructed_momentum_usable = ak.sum(reconstructed_3p[reconstructable_position], axis = 1).to_numpy()
        #total_reconstructed_momentum_usable = total_reconstructed_momentum_all_events[reconstructable_position]
    #np.linalg.norm(total_reconstructed_momentum_all_events, axis = 1) > 0
    
    shortest_distance = get_nearest_approach(data['position'][reconstructable_position], total_reconstructed_momentum_usable)
    
    return shortest_distance, reconstructable_position

def make_plot(distances, weights):
    labels = ['RHN', 'HXX']
    
    fig, ax = plt.subplots(1,1, figsize = (6,4))
    for distance, label, weight in zip(distances, labels, weights):
        print(len(distance),len(weight))
        ax.hist(distance, bins = np.arange(-1, 200, 10), weights = weight, histtype = 'step', label = label)
        
    ax.grid()
    ax.legend(title = 'Model')
    ax.set_xlabel('Distance From IP (m)')
    ax.set_ylabel('Counts')
    ax.set_title('Comparing Closest Approach Distances of Different Models')
    plt.savefig('Closest_Approach.pdf')
    

def main(sim_packs):
    '''
    The main function for the script.  It makes a histogram for the shortest distance
    of reconstructable tracks

    Parameters
    ----------
    sim_packs : list
        Simulation parameters

    Returns
    -------
    distances : list
        Calculated distances from all runs

    '''
    
    distances = []
    
    for i,pack in enumerate(sim_packs):
        
        sim_info = organize_sim_info(*pack)
        vertices =  ad.time_process(ad.do_sim, sim_info)
        
        data = eds.to_data_structure(vertices)
        eds.data_to_pickle(data, f'Simulation_Data_Pack{i}.pickle')
        
        shortest_distance, reconstructable_bool = calculate_shortest_distance(data)
        distances.append(shortest_distance)
        
    make_plot(distances, [data['weight'][reconstructable_bool.to_numpy()],])
    
    return distances
    
if __name__ == '__main__':
    main([['../../SimulationData/ForMixingPlot/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv',
          '../../SimulationData/RHNctauUe.dat', '../../SimulationData/vN_to_all_0.314228_hadron_style.txt'],])