# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 11:08:29 2024

@author: baitk
"""
import sys

import numpy as np
import awkward as ak
import pickle

sys.path.insert(0,'../../FastSim_Additions/')
from Run_Combiner import main as combine 

import matplotlib.pyplot as plt
import boost_histogram as bh
plt.style.use('science')
plt.rcParams.update({'axes.grid':True, 'legend.frameon':True})

from Visibility_Plots import convert_wild_to_label

def initiate_plot():
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6,5))
    
    ax.set_xlabel('Momentum Carried Away From Photons [GeV]')
    ax.set_ylabel('Fraction of Vertices')
    ax.set_yscale('log')
    #ax[1].set_xlabel('Fraction of Photons')
    #x[1].set_ylabel('Fraction of Vertices')
    
    ax.set_xlim(0,120)
    #ax[1].set_xlim(-0.5,18)
    
    return fig, ax

def add_curve(data_column, weights, ax):
    
    photons_involved = data_column > 0
    
    nbins = int((np.max(data_column) - np.min(data_column)) / 5)

    h = bh.Histogram(bh.axis.Regular(nbins,0, np.max(data_column)), storage = bh.storage.Weight())
    h.fill(data_column[photons_involved], weight = weights[photons_involved])

    
    integral = np.sum(h.view().value) * np.diff(h.axes[0].edges)[0]
    
    l = ax.step(h.axes[0].edges[:-1], h.view().value / integral, where = 'post')
    ax.step(h.axes[0].edges[1:], h.view().value / integral, where = 'pre', color = l[0].get_color())
    err = ax.errorbar(h.axes[0].centers, h.view().value / integral, yerr = np.sqrt(h.view().variance)/integral, fmt = '.', color = l[0].get_color())
    
    return err

def get_photon_momenta(data):
    photons = data['daughters'][data['daughters'][:,:,4] == 22]
    
    three_p_from_photons = ak.sum(photons[:,:,1:4], axis = 1)
    square_norm = np.array(ak.sum(three_p_from_photons**2, axis = 1))
    
    return np.sqrt(square_norm)

def tot_events_with_photons(data):
    photons = ak.sum(data['daughters'][:,:,4] == 22,axis = 1)
    
    return np.sum(np.array(photons) != 0) / np.array(photons).shape[0]
    
def get_all_momenta(data):
    three_p = ak.sum(data['daughters'][:,:,1:4], axis = 1)
    square_norm = np.array(ak.sum(three_p**2, axis = 1))
    
    return np.sqrt(square_norm)

def get_fraction_photons(data):
    
    num_photons = ak.sum(data['daughters'][:,:,5] == 22, axis = 1).to_numpy()
    num_particles = ak.sum(data['daughters'][:,:,0] < np.inf, axis = 1).to_numpy()
    
    return num_photons / num_particles 

def main(event_wildcards):
    plot_objects = []
    
    fig, ax = initiate_plot()
    
    for wildcard in event_wildcards:
        data = combine(wildcard)
        
        print('Fraction Photons =', tot_events_with_photons(data))
        
        photon_momentum = get_photon_momenta(data)
        fractional_momentum = photon_momentum #/ get_all_momenta(data)
        
        #fractional_photons = get_fraction_photons(data)
        
        mom_not_nan = ~np.isnan(fractional_momentum)
        #n_not_nan = ~np.isnan(fractional_photons) 
        
        p_frac = add_curve(fractional_momentum[mom_not_nan], data['weight'][mom_not_nan], ax)
        #n_frac = add_curve(fractional_photons[n_not_nan], data['weight'][n_not_nan], ax[1])
        
        plot_objects.append(p_frac)
        
    with open('../LargeRHNSimScripts/sim_SMS_Bmeson_3.06185_10e-10.pickle', 'rb') as f:
        data = pickle.load(f)['Data']
       
    print('Fraction Photons =', tot_events_with_photons(data))
    photon_momentum = get_photon_momenta(data)
    print(photon_momentum)
    fractional_momentum = photon_momentum #/ get_all_momenta(data)
    
    #fractional_photons = get_fraction_photons(data)
    
    mom_not_nan = ~np.isnan(fractional_momentum)
    #n_not_nan = ~np.isnan(fractional_photons) 
    p_frac = add_curve(fractional_momentum[mom_not_nan], data['weight'][mom_not_nan], ax)
    #n_frac = add_curve(fractional_photons[n_not_nan], data['weight'], ax[1])
    
    plot_objects.append(p_frac)
    
    event_wildcards.append('../LargeRHNSimScripts/sim_SMS_Bmeson_3.06185_10e-10.pickle')
    ax.legend([l[0] for l in plot_objects], [convert_wild_to_label(wildcard) for wildcard in event_wildcards])
    #ax[1].legend([l[1] for l in plot_objects], [convert_wild_to_label(wildcard) for wildcard in event_wildcards])
    
    plt.savefig('NumberOfPhotons.pdf')

if __name__ == '__main__':
    runs_to_plot = ['../LargeRHNSimScripts/Finished_Sim_Ue/sim_Ue_*_1.1283_0.00043939705607607863.pickle',
                    '../LargeRHNSimScripts/Finished_Sim_Ue/sim_Ue_*_2.06318_2.275845926074791e-10.pickle',
                    '../LargeRHNSimScripts/Finished_Sim_Ue/sim_Ue_*_3.77268_2.275845926074791e-10.pickle']
     
    #runs_to_plot = []
    with open('../LargeRHNSimScripts/sim_SMS_Bmeson_3.06185_10e-10.pickle', 'rb') as f:
        data = pickle.load(f)['Data']
    main(runs_to_plot)

