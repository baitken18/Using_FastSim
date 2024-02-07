# -*- coding: utf-8 -*-
"""
Creates a plot for number of visible particles per vertex

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

def convert_wild_to_label(wildcard):
    params = wildcard[:-7].split('_')[-4:]
    if params[0] == 'SMS':
        return fr'SMS (m = {params[2]}GeV, $\sin^2\theta$ = {params[3]})'
    else:
        return fr'RHN (m = {params[2]}GeV, $|U|^2$ = {float(params[3]):.4g})'

def initiate_plot():
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (12,5))
    
    ax[0].set_xlabel('Number of Visible Particles')
    ax[0].set_ylabel('Fraction of Vertices')
    ax[1].set_xlabel('Number of Invisible Particles')
    ax[1].set_ylabel('Fraction of Vertices')
    
    ax[0].set_xlim(-0.5,16.5)
    ax[1].set_xlim(-0.5,18.5)
    
    return fig, ax

def add_curve(data_column, weights, ax):
    
    print(data_column)
    
    nbins = np.max(data_column) + 1
        
    print(nbins)
        

    h = bh.Histogram(bh.axis.Regular(nbins,-0.5, nbins - 0.5), storage = bh.storage.Weight())
    h.fill(data_column.to_numpy(), weight = weights)
    
    integral = np.sum(h.view().value)
    
    l = ax.step(h.axes[0].edges, np.append(h.view().value / integral, np.array([h.view().value[-1]/integral])), where = 'post')
    err = ax.errorbar(h.axes[0].centers, h.view().value / integral, yerr = np.sqrt(h.view().variance)/integral, fmt = '.', color = l[0].get_color())
    
    ax.set_xticks(np.arange(0, np.max(data_column) + 1, 2))
    
    return err

def count_visible(data):
    vis_bool = data['daughters'][:,:,-1]
    return ak.sum(vis_bool == 1, axis = 1)

def count_invisible(data):
    vis_bool = data['daughters'][:,:,-1]
    return ak.sum(vis_bool == -1, axis = 1)

def main(event_wildcards):
    plot_objects = []
    
    fig, ax = initiate_plot()
    
    for wildcard in event_wildcards:
        data = combine(wildcard)
        visible = count_visible(data)
        invisible = count_invisible(data)
        
        vis = add_curve(visible, data['weight'], ax[0])
        invis = add_curve(invisible, data['weight'], ax[1])
        
        plot_objects.append((vis,invis))
        
    with open('../LargeRHNSimScripts/sim_SMS_Bmeson_3.06185_10e-10.pickle', 'rb') as f:
        data = pickle.load(f)['Data']
        
    visible = count_visible(data)
    invisible = count_invisible(data)
    
    vis = add_curve(visible, data['weight'], ax[0])
    invis = add_curve(invisible, data['weight'], ax[1])
    
    plot_objects.append((vis,invis))
    
    event_wildcards.append('../LargeRHNSimScripts/sim_SMS_Bmeson_3.06185_10e-10.pickle')
    ax[0].legend([l[0] for l in plot_objects], [convert_wild_to_label(wildcard) for wildcard in event_wildcards])
    ax[1].legend([l[1] for l in plot_objects], [convert_wild_to_label(wildcard) for wildcard in event_wildcards])
    
    plt.savefig('NumberOfVertices_Ver2.pdf')

if __name__ == '__main__':
    runs_to_plot = ['../LargeRHNSimScripts/Finished_Sim_Ue/sim_Ue_*_1.1283_0.00043939705607607863.pickle',
                    '../LargeRHNSimScripts/Finished_Sim_Ue/sim_Ue_*_2.06318_2.275845926074791e-10.pickle',
                    '../LargeRHNSimScripts/Finished_Sim_Ue/sim_Ue_*_3.77268_2.275845926074791e-10.pickle']
    
    main(runs_to_plot)
    '''
    with open('../LargeRHNSimScripts/sim_SMS_Bmeson_3.06185_10e-10.pickle', 'rb') as f:
        data = pickle.load(f)['Data']
        
    vis = count_visible(data)
    invis = count_invisible(data)
    
    fig, ax = initiate_plot()
    add_curve(vis, data['weight'], ax[0])
    add_curve(invis, data['weight'], ax[1])
    plt.show()
    '''