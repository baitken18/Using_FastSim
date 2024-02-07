'''
Running this script will make four plots for the input pickle file:
    1. Opening Angle Histogram
    2. Closest Approach Histogram
    3. Scatter plot of closest approach vs opening angle
    4. Opening angle vs 3p magnitude
    
Call with $make_plots sample_data.pickle
'''
import sys

import pickle

import numpy as np
import awkward as ak

import boost_histogram as bh
import matplotlib.pyplot as plt
import matplotlib.colors as colours

plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6,4), 'legend.frameon':True, 'font.size':14})


sys.path.insert(0,'../MomentumExtrapolation/')
from Opening_Angle_Histogram import calculate_opening_angle
from Closest_Approach_Histogram import calculate_shortest_distance

def format_plot(ax, axis_labels, log):
    ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000 fb^{-1}$)')
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    #ax.grid()
    ax.legend()
    if log:
        ax.set_yscale('log')
        
def add_endpoints(boost_hist):
    
    start = boost_hist.axes[0].centers[0] - np.diff(boost_hist.axes[0].centers)[0]
    end = boost_hist.axes[0].centers[-1] + np.diff(boost_hist.axes[0].centers)[0]
    extended_bins = np.array([start, *boost_hist.axes[0].centers, end])
    
    extended_values = np.array([0, *boost_hist.view().value, 0])
    extended_weight = np.array([0, *boost_hist.view().variance, 0])
    
    return extended_bins, extended_values, extended_weight
    
def make_opening_histogram(data, params, outfile = None, ax = None, log = True):
    flavour, mass, mixing_sq = params
    
    angles, weight_bool = calculate_opening_angle(data)
    weights = data['weight'][weight_bool]
    
    nbins = int(np.max(angles) / 0.1 + 1)
    
    h = bh.Histogram(bh.axis.Regular(nbins,0, np.max(angles)), storage = bh.storage.Weight())
    h.fill(angles, weight = weights)
    
    plot_bins, plot_vals, plot_weights = add_endpoints(h)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        add_labs = True
    else:
        add_labs = False
    
    l = ax.step(plot_bins / np.pi, plot_vals, where = 'mid', 
                #label = fr'RHN (Mass = {mass}, $\vert U_\{flavour[1:]} \vert^2$ = {mixing_sq:.3g})')
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.3g})')
    ax.errorbar(plot_bins[1:-1] / np.pi, plot_vals[1:-1], yerr = np.sqrt(plot_weights[1:-1]), fmt = '.', color = l[0].get_color())
    
    if flavour == 'SMS':
        l[0].set_label(fr'SMS (Mass = {mass}, $\sin^2\theta = ${mixing_sq}')
    '''
    l = ax.step(h.axes[0].centers/np.pi, h.view().value, where = 'mid', 
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.5g})')
    ax.errorbar(h.axes[0].centers/np.pi, h.view().value, yerr = np.sqrt(h.view().variance), fmt = '.', color = l[0].get_color())
    '''
    if add_labs:
        format_plot(ax, [r'Maximum Opening Angle [$\pi$ rad]', 'Number of Vertices'], log)
    
    if outfile != None:
        plt.savefig(outfile)
        plt.close()
    
    return ax

def make_closest_histogram(data, params, outfile = None, ax = None, log = True):
    flavour, mass, mixing_sq = params
    
    distances, reconstructable_bool = calculate_shortest_distance(data)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    
    nbins = int((np.max(distances) - np.min(distances)) / 10 + 1) #division is the width of each bin
    
    h = bh.Histogram(bh.axis.Regular(nbins, 0, np.max(distances)),\
                     storage = bh.storage.Weight())
    h.fill(distances, weight = weights)
    
    plot_bins, plot_vals, plot_weights = add_endpoints(h)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Distance of Closest Approach', ['Shortest Distance [m]', 'Number of Vertices'])
        if log:
            ax.set_yscale('log')
    
    l = ax.step(plot_bins, plot_vals, where = 'mid', 
                #label = fr'RHN (Mass = {mass}, $\vert U_\{flavour[1:]} \vert^2$ = {mixing_sq:.3g})')
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.3g})')
    ax.errorbar(plot_bins[1:-1], plot_vals[1:-1], yerr = np.sqrt(plot_weights[1:-1]), fmt = '.', color = l[0].get_color())
    '''
    l = ax.step(h.axes[0].centers, h.view().value, where = 'mid',
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.5g})')
    ax.errorbar(h.axes[0].centers, h.view().value, yerr = np.sqrt(h.view().variance), fmt = '.', color = l[0].get_color())
    '''
    
    if outfile != None:
        plt.savefig(outfile)
        
        plt.close()
        
    return ax

def make_phi_angle_histogram(data, params, outfile = None, ax = None, log = False):
    distances, reconstructable_bool = calculate_shortest_distance(data)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    radii = np.linalg.norm(data['position'][reconstructable_bool.to_numpy()], axis = 1)
    
    flavour, mass, mixing_sq = params
    
    phi = np.arctan(distances / radii)
    
    h = bh.Histogram(bh.axis.Regular(15,0, np.pi/2), storage = bh.storage.Weight())
    h.fill(phi, weight = weights)
    
    plot_bins, plot_vals, plot_weights = add_endpoints(h)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Track Correction Angle', [r'Angle [$\pi$ rad]', 'Number of Vertices'])
        if log:
            ax.set_yscale('log')
            
    l = ax.step(plot_bins / np.pi, plot_vals, where = 'mid', 
                #label = fr'RHN (Mass = {mass}, $\vert U_\{flavour[1:]} \vert^2$ = {mixing_sq:.3g})')
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.3g})')
    ax.errorbar(plot_bins[1:-1] / np.pi, plot_vals[1:-1], yerr = np.sqrt(plot_weights[1:-1]), fmt = '.', color = l[0].get_color())
    '''
    #Change to shift errorbars instead
    l = ax.step(h.axes[0].centers/np.pi, h.view().value, where = 'mid',
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.5g})')
    ax.errorbar(h.axes[0].centers/np.pi, h.view().value, yerr = np.sqrt(h.view().variance), fmt = '.', color = l[0].get_color())
    '''
    
    if outfile != None:
        plt.savefig(outfile)
        plt.close()
    
    return ax
    
    

def compare_opening_magnitude(data, params, outfile = None, ax = None):
    flavour, mass, mixing_sq = params
    
    angles, event_bool = calculate_opening_angle(data)
    mag_3p = np.linalg.norm(data['momentum'][:,1:], axis = 1)[event_bool]
    weights = data['weight'][event_bool]
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000$ fb$^{-1}$)')
        ax.set_xlabel(r'Opening Angle [$\pi$ rad]')
        ax.set_ylabel(r'$\Vert p \Vert$ [GeV]')
        #ax.grid()
        #add_legend = True
    
    angle_bin_max = np.max(angles / np.pi)
    mom_bin_max = 200
    
    h = ax.hist2d(angles / np.pi, mag_3p, bins = [np.linspace(0,angle_bin_max,20), np.linspace(0,mom_bin_max,20)], 
                  weights = weights, label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.5g})',
                  norm = colours.LogNorm(clip = True))
    fig.colorbar(h[3], ax = ax, label = 'Number of Vertices')
    
    #if add_legend:
        #ax.legend()
    
    if outfile != None:
        plt.savefig(outfile)
        plt.show()
        
        plt.close()
    
    return ax
    
def compare_opening_closest(data, params, outfile = None, ax = None):
    flavour, mass, mixing_sq = params
    
    angles, at_least_2 = calculate_opening_angle(data)
    
    for key in data.keys():
        if type(data[key]) != list:
            data[key] = data[key][at_least_2]
    
    distance, _ = calculate_shortest_distance(data)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000$ fb$^{-1}$)')
        ax.set_xlabel(r'Opening Angle [$\pi$ rad]')
        ax.set_ylabel(r'Shortest Distance to CMS IP')
        ax.grid()
        #add_legend = True
    
    h = ax.hist2d(angles / np.pi, distance, bins = [20,20], weights = data['weight'],
                  label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq:.5g})',
                  norm = colours.LogNorm(clip = True))
    fig.colorbar(h[3], ax = ax, label = 'Number of Vertices')
    
    #if add_legend:
        #ax.legend()
    
    if outfile != None:
        plt.savefig(outfile)
        plt.show()
        plt.close()
        
    return ax

def main(filename = None, outfile_prefix = None):
    
    if filename == None:
        filename, outfile_prefix = sys.argv[1:]
        
    mass, mixing = map(float, filename[:-7].split('_')[-2:])
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)['Data']
    
    make_opening_histogram(data, filename, outfile_prefix + 'opening_angle.pdf', log = False)
    make_opening_histogram(data, filename, outfile_prefix + 'opeing_angle_log.pdf')
    make_closest_histogram(data, filename, outfile_prefix + 'closest.pdf', log = False)
    make_closest_histogram(data, filename, outfile_prefix + 'closest_log.pdf')
    compare_opening_closest(data,filename, outfile_prefix + 'opening_closest.pdf')
    compare_opening_magnitude(data, filename, outfile_prefix + 'opening_magnitude.pdf')
    
if __name__ == '__main__':
    main()
    '''
    filename = 'sim_0.1_2.335721469090121e-06.pickle'

    with open(filename, 'rb') as f:
        data = pickle.load(f)['Data']
        
    make_opening_histogram(data)
    make_opening_histogram(data, log = False)
    
    make_closest_histogram(data)
    make_closest_histogram(data, log = False)
    '''
