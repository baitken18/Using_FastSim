import sys

import numpy as np
import awkward as ak

sys.path.insert(0,'../MomentumExtrapolation/')
from Opening_Angle_Histogram import calculate_opening_angle
from Closest_Approach_Histogram import calculate_shortest_distance, get_reconstructable

sys.path.insert(0, '../NeutrinoTools')
from Neutrino_Plots import calculate_reconstruction

def param_weight_order(tup):
    return tup[0]
    
def make_cut(parameter, weight, confidence):
    
    unit_weight = weight / np.sum(weight)
    
    sorted_param_weight = np.array(sorted([np.array((p,w)) for p,w in zip(parameter,unit_weight)],
                                 key = param_weight_order))
    
    greater_mask = np.cumsum(sorted_param_weight.T[1]) >= confidence
    
    return sorted_param_weight[greater_mask][0,0] 

def check_cut(parameter_values, cut, confidence):
    
    unit_weights = parameter_values[0] / np.sum(parameter_values[0])
    
    opening_bool = parameter_values[1] < cut[0]
    recon_bool = parameter_values[2] < cut[1]
    shortest_bool = parameter_values[3] < cut[2]
    
    stage1 = opening_bool & recon_bool
    stage2 = stage1 & shortest_bool
    
    perc_events = np.sum(unit_weights[stage2])
    
    return perc_events > confidence #and perc_events < confidence + 0.003
    

def fetch_values(data, track_thresh):
    
    opening,event_bool = calculate_opening_angle(data, track_thresh)
    recon = calculate_reconstruction(data, track_thresh)[0]
    shortest = calculate_shortest_distance(data, track_thresh)[0]
    
    return np.array([data['weight'][event_bool], opening, recon, shortest])

def get_ranges(parameter_values, prec, con):
    
    min_vals = np.array([make_cut(p,parameter_values[0], confidence = con) 
                         for p in parameter_values[1:]])
    max_vals = np.max(parameter_values[1:], axis = 1)
    
    return [np.arange(mi, ma, p) for mi,ma,p in zip(min_vals, max_vals, prec)]

def run_algorithm(parameter_values, track_thresh = 2, prec = (0.01, 0.01, 1), con = 0.95):
    
    opening_range, recon_range, short_range = get_ranges(parameter_values, prec, con)
    
    potential_cuts = []
    
    for open_cut in opening_range:
        for recon_cut in recon_range:
            for short_cut in short_range:
                packet = (open_cut, recon_cut, short_cut)
                if check_cut(parameter_values, packet, con):
                    potential_cuts.append(packet)
                    
    return np.array(potential_cuts)
    
def main(data, track_thresh = 2, prec = (0.1, 0.1,5), con = 0.95):
    
    parameters = fetch_values(data, track_thresh)
    
    return run_algorithm(parameters, track_thresh, prec, con)
    