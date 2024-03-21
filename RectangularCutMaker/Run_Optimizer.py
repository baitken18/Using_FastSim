import sys
from CutAlgorithm import main as get_cuts
import pickle 
import glob
import boost_histogram as bh

import numpy as np
import awkward as ak

sys.path.insert(0, '../../FastSim_Additions/')
from Run_Combiner import main as combine

sys.path.insert(0,'../MomentumExtrapolation/')
from Opening_Angle_Histogram import calculate_opening_angle
from Closest_Approach_Histogram import calculate_shortest_distance

sys.path.insert(0, '../NeutrinoTools')
from Neutrino_Plots import calculate_reconstruction, make_1d_histogram

sys.path.insert(0, '../FinalPlots')
from FinalPlotGenerator import collect_Bgd_runs
'''
def collapse(data_dict):

    all_data = data_dict[0]
    for run in range(1,len(data_dict.keys())):
        print(f'Collapsing Run {run}')
        for quantity in data_dict[run].keys():
            if quantity == 'daughters':
                all_data[quantity] = ak.flatten(ak.Array([all_data[quantity],
                                                         data_dict[run][quantity]]))
                
            else:
                all_data[quantity] = np.append(all_data[quantity],
                                               data_dict[run][quantity])

    return all_data
'''

def get_quantities(data, track_thresh):
    opening, weight_bool = calculate_opening_angle(data, track_thresh)
    recon = calculate_reconstruction(data, track_thresh)[0]
    shortest = calculate_shortest_distance(data, track_thresh)[0]
    return data['weight'][weight_bool], opening, recon, shortest
                    
def make_bgd_hists(filenames, track_thresh = 2):
    
    h = bh.Histogram(bh.axis.Regular(100, 0, np.pi), bh.axis.Regular(100,0,np.pi),
                     bh.axis.Regular(100, 0, 200), storage = bh.storage.Weight())
    
    for fname in filenames:
        print(f'Starting file -> {fname}')
        with open(fname, 'rb') as f:
            data = pickle.load(f)
            
        quantities = get_quantities(data, track_thresh)
        print('Quantities Calculated')
        h.fill(*quantities[1:], weight = quantities[0])
        print('Histogram Filled')
        
    with open('BgdHistDV3.pickle', 'wb') as f:
        pickle.dump(h,f)
            
    return h
        

def bgd_cut(data, cut, track_thresh):
    
    opening = np.where(data.axes[0].edges < cut[0])[0][-1]
    recon = np.where(data.axes[1].edges < cut[1])[0][-1]
    short = np.where(data.axes[2].edges < cut[2])[0][-1]
    
    sum_events = np.sum(data.view().variance[:opening,:recon,:short])
    tot_events = np.sum(data.view().variance)
    
    return 1 - (sum_events / tot_events)

def get_optimal(RHN_data, Bgd_data, track_thresh = 2, con = 0.95):
    print(f'Started {con}')
    
    potential_cuts = get_cuts(RHN_data, track_thresh = track_thresh, prec = (0.1, 0.1,5), con = con)
    
    if len(potential_cuts) == 0:
        return np.zeros(3), 0
    
    optimal_cut = potential_cuts[0]
    optimal_perc = bgd_cut(Bgd_data, potential_cuts[0], track_thresh)
    
    for cut in potential_cuts[1:]:
        background_perc = bgd_cut(Bgd_data, cut, track_thresh)
        
        if background_perc > optimal_perc:
            optimal_cut = cut
            optimal_perc = background_perc
    
    return optimal_cut, optimal_perc

def main(RHN_params, track_thresh = 2, con = 0.95):
    RHN_data = combine(RHN_params)
    
    #Bgd_data = make_bgd_hists(glob.glob('../NeutrinoResults/*.pickle'))
    with open('BgdHist.pickle', 'rb') as f: 
        Bgd_data = pickle.load(f)
        print('Data Loaded')
    
    return get_optimal(RHN_data, Bgd_data, track_thresh, con) 

if __name__ == '__main__':
    #filenames = glob.glob('../NeutrinoResults/*.pickle')
    #hists = make_bgd_hists(filenames)
    RHN_params = ('Ue',4.08883,5.1794746792312124e-08)
    
    optimal_cut, optimal_perc = main(RHN_params, con = 0.95)
    print(optimal_cut)
