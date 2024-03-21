# -*- coding: utf-8 -*-
"""
The module containing the foundation of all the neutrino plots
"""
import sys

import numpy as np
import awkward as ak

import boost_histogram as bh
import matplotlib.pyplot as plt
import matplotlib.colors as colours

#plt.style.use('science')
#plt.rcParams.update({'figure.figsize':(6,4), 'legend.frameon':True, 'font.size':14})

sys.path.insert(0,'../MomentumExtrapolation/')
from Opening_Angle_Histogram import calculate_opening_angle
from Closest_Approach_Histogram import calculate_shortest_distance, get_reconstructable

def make_label(params):
    if len(params) == 3:
        return fr'Background $\nu_\mu$ ({params[1]}-{params[2]}GeV)'
        #return fr'RHN (M={params[1]}GeV, $\vert U_e\vert^2$={params[2]:.3g})'
    else: 
        return 'Total Neutrino Background'

def make_1d_histogram(data_col, data_weight, bin_params):
    
    h = bh.Histogram(bh.axis.Regular(*bin_params), storage = bh.storage.Weight())
    h.fill(data_col, weight = data_weight)
    
    return h

def draw_1d_histogram(h, ax, label, x_scaling = 1, normalize = True):
    
    if normalize:
        norm = np.sum(h.view().value)
        print('Number of Events', norm)
        h.view().value = h.view().value / norm
        h.view().variance = h.view().variance / norm**2
    
    l = ax.step(h.axes[0].edges / x_scaling, np.append(np.zeros(1), h.view().value),
                where = 'pre', label = label)
    ax.step(h.axes[0].edges / x_scaling, np.append(h.view().value, np.zeros(1)),
            where = 'post', color = l[0].get_color())
    ax.errorbar(h.axes[0].centers / x_scaling, h.view().value, yerr = np.sqrt(h.view().variance), 
                fmt = 'none', color = l[0].get_color())

def add_opening_curve(data, params, ax):
    
    angles, weight_bool = calculate_opening_angle(data, 3)
    weights = data['weight'][weight_bool]
    
    nbins = 20
    
    hist = make_1d_histogram(angles, weights, (nbins, 0, np.pi))
    draw_1d_histogram(hist,ax, make_label(params), x_scaling = np.pi)
    
    return ax

def cumulative_opening_curve(data_sets, ax):
    
    val_storage = np.array([])
    weight_storage = np.array([])
    
    for key in data_sets.keys():
        
        angles, weight_bool = calculate_opening_angle(data_sets[key], 2)
        weights = data_sets[key]['weight'][weight_bool]
        
        val_storage = np.append(val_storage, angles)
        weight_storage = np.append(weight_storage, weights)
    
    nbins = 20
    
    hist = make_1d_histogram(val_storage, weight_storage, (nbins, 0, np.pi))
    draw_1d_histogram(hist,ax, make_label(''), x_scaling = np.pi)
    
    return ax

def add_closest_curve(data, params, ax):
    
    distances, reconstructable_bool = calculate_shortest_distance(data,2)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    
    nbins = 40
    hist = make_1d_histogram(distances, weights, (nbins, 0, 200))
    draw_1d_histogram(hist, ax, make_label(params))
    
    return ax 

def cumulative_closest_curve(data_sets, ax):
    
    val_storage = np.array([])
    weight_storage = np.array([])
    
    for key in data_sets.keys():
        
        distances, weight_bool = calculate_shortest_distance(data_sets[key],2)
        weights = data_sets[key]['weight'][weight_bool.to_numpy()]
        
        val_storage = np.append(val_storage, distances)
        weight_storage = np.append(weight_storage, weights)
    
    nbins = 40
    hist = make_1d_histogram(val_storage, weight_storage, (nbins, 0, 200))
    draw_1d_histogram(hist, ax, make_label(''))
    
    return ax

def get_total_3p(data, reconstructable):
    reconstructed_3p = ak.fill_none(ak.mask(data['daughters'][:,:,1:4], reconstructable),
                                    ak.Array([0,0,0]),axis = 1)
    three_p = ak.sum(reconstructed_3p, axis = 1)
    
    return three_p

def calculate_reconstruction(data, track_thresh = 2):
    reconstructable = get_reconstructable(data)
    recon_events = ak.sum(reconstructable, axis = 1) >= track_thresh
    total_3p = get_total_3p(data, reconstructable)[recon_events]
    try:
        positions = data['position'][recon_events]
    except IndexError:
        positions = data['position'].reshape(100000, 3)[recon_events]
    
    dot_p = np.sum(total_3p * positions, axis = 1)
    cos_arg = dot_p / (np.linalg.norm(total_3p, axis = 1) * np.linalg.norm(positions, axis = 1))
    cos_arg = np.array(cos_arg)
    cos_arg[np.where(cos_arg > 1)[0]] = 1
    phi = np.arccos(cos_arg)
    
    
    return phi, recon_events

def add_reconstruction_curve(data, params, ax):
    '''
    distances, reconstructable_bool = calculate_shortest_distance(data)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    radii = np.linalg.norm(data['position'][reconstructable_bool.to_numpy()], axis = 1)
    phi = np.arcsin(distances / radii)
    '''
    phi,reconstructable = calculate_reconstruction(data)
    weights = data['weight'][reconstructable]
    
    nbins = 40
    hist = make_1d_histogram(phi, weights, (nbins, np.min(phi), np.max(phi)))
    draw_1d_histogram(hist, ax, make_label(params), x_scaling = np.pi)
    
    return ax 

def cumulative_reconstruction_curve(data_sets, ax):
    
    val_storage = np.array([])
    weight_storage = np.array([])
    
    for key in data_sets.keys():
        '''
        distances, weight_bool = calculate_shortest_distance(data_sets[key])
        weights = data_sets[key]['weight'][weight_bool.to_numpy()]
        radii = np.linalg.norm(data_sets[key]['position'][weight_bool.to_numpy()], axis = 1)
        '''
        phi,reconstructable = calculate_reconstruction(data_sets[key],2)
        weights = data_sets[key]['weight'][reconstructable]
        
        val_storage = np.append(val_storage, phi)
        weight_storage = np.append(weight_storage, weights)
    
    nbins = 40
    hist = make_1d_histogram(val_storage, weight_storage, (nbins, 0, np.pi))
    draw_1d_histogram(hist, ax, make_label(''), x_scaling = np.pi)
    
    return ax

def compare_opening_closest(data, params, fig, ax):
    
    angles, at_least_2 = calculate_opening_angle(data)
    for key in data.keys():
        if type(data[key]) != list:
            data[key] = data[key][at_least_2]
    distance, _ = calculate_shortest_distance(data)
    
    h = ax.hist2d(angles / np.pi, distance, bins = [20,20], weights = data['weight'],
                  norm = colours.LogNorm(clip = True))
    fig.colorbar(h[3], ax = ax, label = 'Number of Vertices')
    
    return ax
    

    
    
    
    
    
    
    