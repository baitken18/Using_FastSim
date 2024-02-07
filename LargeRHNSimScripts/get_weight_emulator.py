# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 20:42:59 2024

@author: baitk
"""

import numpy as np

def get_weight_emulate(decay_vol_hits):
    
    #decay_vol_hits.sort(key = np.linalg.norm)
    
    L1, L2 = np.linalg.norm(decay_vol_hits, axis = 1)
    
    const = 800
    
    unif = np.random.uniform(1 - np.exp(-L1/const), 1 - np.exp(-L2/const))
    exp = - const * np.log(1 - unif)
    
    decay_position = np.array(decay_vol_hits[0])/L1 * exp
    
    return decay_position

def get_trajectories(n):
    angles = np.linspace(np.pi/6, np.pi / 3, n)
    
    tan_angle = np.tan(angles)
    
    all_decay_vol_hits = np.zeros((n,2,2))
    
    all_decay_vol_hits[angles < np.pi/4] = np.array([np.array([[50 / ang, 50],[50*np.sqrt(3), 50 * np.sqrt(3) * ang]]) 
                                                      for ang in tan_angle[angles < np.pi / 4]])
    
    all_decay_vol_hits[angles > np.pi/4] = np.array([np.array([[50, 50/ang],[50*np.sqrt(3)*ang, 50 * np.sqrt(3)]]) 
                                                      for ang in tan_angle[angles < np.pi / 4]])
    
    return all_decay_vol_hits

def main(n = 30):
    all_decay_vol_hits = get_trajectories(n)
    
    positions = [get_weight_emulate(hit) for hit in all_decay_vol_hits]
    
    Ls = np.linalg.norm(all_decay_vol_hits, axis = 2)
    
    return np.array(positions), Ls
    
    
                                                      
                                                      
    
    