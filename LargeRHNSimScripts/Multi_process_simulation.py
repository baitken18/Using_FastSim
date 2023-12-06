# -*- coding: utf-8 -*-
"""
Script will be used to run large simulations for RHN process
"""
import sys
import os
from multiprocessing import Pool

import numpy as np
import timeit

sys.path.insert(0, '../../FastSim_Additions/')
import Additions as ad
import Event_Data_Structure as eds

sys.path.insert(0, '../../MATHUSLA_FastSim/')
import DetectorSimulation.llp_gun_new as lg
import Helpers.functions as hf

class sim_pack:
    '''Sim_pack is used to parameterize each simulation'''
    def __init__(self, fv_file:str, prod_file:str, length_file:str, mixing_range:tuple, sim_size:int):
        self.fv_path = os.path.join(os.getcwd(), fv_file)
        self.prod_path = os.path.join(os.getcwd(), prod_file)
        self.length_path = os.path.join(os.getcwd(), length_file)
        self.mass = float(fv_file.split('_')[-1][:-4])
        self.mixing_range = (int(mixing_range[0]), int(mixing_range[1]))
        self.mixing_res = 20
        self.sim_size = sim_size
    
    def get_mixing_sq(self):
        '''Takes a log scale for mixing_sq based on given range'''
        return np.logspace(self.mixing_range[0],self.mixing_range[1], self.mixing_res)
    
    def get_proper_lifetime(self, mixing_sq):
        '''Gives the proper lifetimes based off mass and range of mixing angles'''
        ctaus = ad.get_ctaus(mass = self.mass, mixing_sq = mixing_sq, length_file = self.length_path)
        return ctaus
    
    def get_outfile(self, mixing_sq):
        '''Outfile name for a dataset based on parameter space'''
        return f'sim_{self.mass}_{mixing_sq}.pickle'

def make_sim_packs(param_file):
    '''Creates list of sim_pack objects based off of the param_file passed
    The param file will be made by param_file_writer'''
    sim_packs = []
    param_content = np.loadtxt(param_file, delimiter = ' ', dtype = str)
    for param_set in param_content:
        print(param_set)
        param_set = list(param_set)
        param_set[3] = tuple(param_set[3][1:-1].split(','))
        param_set[4] = int(param_set[4])
        sim_packs.append(sim_pack(*param_set))
    
    return sim_packs 
    
def single_do_sim(sim_pack, mixing_sq):
    print(f'Starting Mass: {sim_pack.mass} | Mixing_sq: {mixing_sq}')

    #Sets up the detector object and gets relevant angles
    detector_benchmark = ad.initiate_detector()
    phi_min, phi_max, theta_min, theta_max = hf.get_detector_angles(detector_benchmark)
    
    storage = []
    ctau = sim_pack.get_proper_lifetime(mixing_sq)
        
    #Takes the four vectors from events and reweights points
    all_llp_4p = np.array(hf.read_vectors(sim_pack.fv_path, sim_pack.sim_size))[:,1:]
    all_llp_weights = np.array(hf.read_vectors(sim_pack.fv_path, sim_pack.sim_size))[:,0]
    all_llp_weights = ad.rescale_weights(all_llp_weights, phi_min, phi_max, mixing_sq, sim_pack.fv_path)
    
    #Loops through four vectors and their weights
    for llp_4p,w in zip(all_llp_4p,all_llp_weights):
        
        llp_theta = hf.get_theta(llp_4p[1:])
        
        #Checks if particle theta in detector
        if (llp_theta < theta_max) and (llp_theta > theta_min):
            detector_benchmark.clear_detector()
            
            rotated_four_p = hf.deal_with_phi(llp_4p, phi_min, phi_max)
            
            pack = hf.get_weight(rotated_four_p, sim_pack.mass, ctau, detector_benchmark)
            
            if pack is not None:
                p_decay, p_pos, boost = pack
                
                #This step gets the vertex object
                llp_vertex = lg.get_llp('hadronic', sim_pack.mass, p_pos, boost, sim_pack.prod_path)
            
                if llp_vertex is not None:
                    
                    detector_benchmark.new_vertex_event(llp_vertex)
                    
                    storage.append((w,detector_benchmark.return_current_vertex()))
                        
    return storage


def multi_do_sim(sim_pack):

    start = timeit.default_timer()
    
    #Loops over mixing angles from the sim_pack
    for mixing_sq in sim_pack.get_mixing_sq():
        
        #Does a single simulations for a given sim_pack and mixing_sq
        vertices = single_do_sim(sim_pack, mixing_sq)
        
        #Exports the data into a pickle file
        data = eds.to_data_structure(vertices)
        eds.data_to_pickle(data, sim_pack.get_outfile(mixing_sq))
        
    end = timeit.default_timer()
    print(f'Mass: {sim_pack.mass}GeV | Time: {end - start}')
        

def main(param_file = None, cores = 4):
    '''
    Runs entire script to do multiprocess of events in the param_file passed to it

    Parameters
    ----------
    param_file : str, optional
        Parameter file made by param_file_writer The default is None.
    cores : int, optional
        Number of cores to run the simulation. The default is 4.

    Returns
    -------
    None.

    '''
    if param_file == None:
        _, cores, param_file = sys.argv
        
    sim_packs = make_sim_packs(param_file)
    
    p = Pool(int(cores))
    
    p.map(multi_do_sim, sim_packs)

if __name__ == '__main__':
    main()

