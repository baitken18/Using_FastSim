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

def organize_sim_info(fv_file, length_file, product_file, param_card_file = '../../MATHUSLA_FastSim/param_card_CDR.txt'):
    #Relative Paths to important Information
    fv_path = os.path.join(os.getcwd(), fv_file)
    len_path = os.path.join(os.getcwd(), length_file)
    prod_path = os.path.join(os.getcwd(), product_file)

    mass = float(fv_file.split('_')[-1][:-4])
    mixing = 0.004
    ctau = ad.get_ctaus(mass = mass, mixing = mixing, length_file = len_path)
    detector_benchmark = ad.initiate_detector('../../MATHUSLA_FastSim/param_card_CDR.txt')
    
    #The last integer is the simulation size
    return [fv_path, prod_path, ctau, mass, detector_benchmark, 20000]

def get_reconstructable(data, hit_threshold = 4):
    enough_hits = data['daughters'][:,:,5] >= hit_threshold
    not_invisible = data['daughters'][:,:,6] > -1
    
    return enough_hits & not_invisible

def separate_dot(x,y):
    return np.sum(x * y, axis = 1)

def proj(x,y):
    return separate_dot(x,y) * y.T / np.linalg.norm(y, axis = 1)**2

def get_nearest_approach(position, momentum):
    projection = position - proj(position, momentum).T
    return np.linalg.norm(projection, axis = 1)

def calculate_shortest_distance(data):
    
    #These are boolean masks which tell when events are reconstructable
    #Second line checks that the vertex has at least one reconstructed daughter
    reconstructable_daughters = get_reconstructable(data)
    reconstructable_position = ak.sum(reconstructable_daughters, axis = 1) >= 1
    
    reconstructed_3p = ak.fill_none(ak.mask(data['daughters'][:,:,1:4], reconstructable_daughters),
                                    ak.Array([0,0,0]),axis = 1)
    
    total_reconstructed_momentum_all_events = ak.sum(reconstructed_3p, axis = 1).to_numpy()
    total_reconstructed_momentum_usable = total_reconstructed_momentum_all_events[np.linalg.norm(total_reconstructed_momentum_all_events, axis = 1) > 0]
    
    shortest_distance = get_nearest_approach(data['position'][reconstructable_position], total_reconstructed_momentum_usable)
    
    return shortest_distance

def make_plot(distances):
    labels = ['RHN', 'HXX']
    
    fig, ax = plt.subplots(1,1, figsize = (6,4))
    for distance, label in zip(distances, labels):
        ax.hist(distance, bins = np.arange(-1, 200, 10), histtype = 'step', label = label)
    
    ax.grid()
    ax.legend(title = 'Model')
    ax.set_xlabel('Distance From IP')
    ax.set_ylabel('Counts')
    ax.set_title('Comparing Closest Approach Distances of Different Models')
    plt.savefig('Closest_Approach.pdf')
    

def main(sim_packs):
    
    distances = []
    
    for i,pack in enumerate(sim_packs):
        
        sim_info = organize_sim_info(*pack)
        vertices =  ad.time_process(ad.do_sim, sim_info)
        
        data = eds.to_data_structure(vertices)
        eds.data_to_pickle(data, f'Simulation_Data_Pack{i}.pickle')
        
        shortest_distance = calculate_shortest_distance(data)
        distances.append(shortest_distance)
        
    make_plot(distances)
    
if __name__ == '__main__':
    main([['../../SimulationData/ForMixingPlot/RHN_Ue_LLPweight4vectorBmesonlist_mN_0.316228.csv',
          '../../SimulationData/RHNctauUe.dat', '../../SimulationData/vN_to_all_0.314228_hadron_style.txt'],])