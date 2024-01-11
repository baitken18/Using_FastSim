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
plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6,4), 'legend.frameon':True, 'font.size':14})


sys.path.insert(0,'../MomentumExtrapolation/')
from Opening_Angle_Histogram import calculate_opening_angle
from Closest_Approach_Histogram import calculate_shortest_distance

def format_plot(ax, axis_labels, log):
    ax.set_title(r'MATHUSLA100 No Walls Geometry (Luminosity = 3000 fb$^{-1})')
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.grid()
    ax.legend()
    if log:
        ax.set_yscale('log')
    
def make_opening_histogram(data, filename, outfile = None, ax = None, log = True):
    flavour, _, mass, mixing_sq = filename[:-7].split('_')[-4:]
    
    angles, weight_bool = calculate_opening_angle(data)
    weights = data['weight'][weight_bool]
    
    nbins = int(np.max(angles) / 0.1 + 1)
    
    h = bh.Histogram(bh.axis.Regular(nbins,0, np.max(angles)), storage = bh.storage.Weight())
    h.fill(angles, weight = weights)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
    
    #Change to shift errorbars instead
    l = ax.step(h.axes[0].centers/np.pi, h.view().value, where = 'mid', 
                label = fr'RHN (Mass = {mass}, $\vert {flavour} \vert^2$ = {mixing_sq})')
    ax.errorbar(h.axes[0].centers/np.pi, h.view().value, yerr = np.sqrt(h.view().variance), fmt = '.', color = l[0].get_color())
    
    if ax == None:
        format_plot(ax, [r'Maximum Opening Angle [$\pi$ rad]', 'Number of Vertices'], log)
    
    if outfile != None:
        plt.savefig(outfile)
        plt.close()
    
    return ax

def make_closest_histogram(data, outfile = None, ax = None, log = True):
    distances, reconstructable_bool = calculate_shortest_distance(data)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    
    nbins = int((np.max(distances) - np.min(distances)) / 10 + 1) #division is the width of each bin
    
    h = bh.Histogram(bh.axis.Regular(nbins, 0, np.max(distances)),\
                     storage = bh.storage.Weight())
    h.fill(distances, weight = weights)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Distance of Closest Approach', ['Shortest Distance [m]', 'Number of Vertices'])
        if log:
            ax.set_yscale('log')
    
    l = ax.step(h.axes[0].centers, h.view().value, where = 'mid')
    ax.errorbar(h.axes[0].centers, h.view().value, yerr = np.sqrt(h.view().variance), fmt = '.', color = l[0].get_color())
    
    if outfile != None:
        plt.savefig(outfile)
        
        plt.close()
        
    return ax

def make_phi_angle_histogram(data, outfile = None, ax = None, log = False):
    distances, reconstructable_bool = calculate_shortest_distance(data)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    radii = data['position'][reconstructable_bool.to_numpy()]
    
    phi = np.arctan(distances / radii)
    
    h = bh.Histogram(bh.axis.Regular(30,-np.pi, np.pi), storage = bh.storage.Weight())
    h.fill(phi, weight = weights)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Track Correction Angle', [r'Angle [$\pi$ rad]', 'Number of Vertices'])
        if log:
            ax.set_yscale('log')
    
    #Change to shift errorbars instead
    l = ax.step(h.axes[0].centers/np.pi, h.view().value, where = 'mid')
    ax.errorbar(h.axes[0].centers/np.pi, h.view().value, yerr = np.sqrt(h.view().variance), fmt = '.', color = l[0].get_color())
    
    if outfile != None:
        plt.savefig(outfile)
        plt.close()
    
    return ax
    
    

def compare_opening_magnitude(data, outfile = None, ax = None):
    angles, event_bool = calculate_opening_angle(data)
    mag_3p = np.linalg.norm(data['momentum'][:,1:], axis = 1)[event_bool]
    weights = data['weight'][event_bool]
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Comparing 3-Momentum to Opening Angle', [r'Opening Angle [$\pi$ rad]', r'$\Vert p \Vert$ [GeV]'])
        ax.grid()
    
    h = ax.hist2d(angles / np.pi, mag_3p, bins = [20,20], weights = weights)
    fig.colorbar(h[3], ax = ax, label = 'Number of Vertices')
    
    if outfile != None:
        plt.savefig(outfile)
        
        plt.close()
    
    return ax
    
def compare_opening_closest(data, outfile = None, ax = None):
    angles, at_least_2 = calculate_opening_angle(data)
    
    for key in data.keys():
        if type(data[key]) != list:
            data[key] = data[key][at_least_2]
    
    distance, _ = calculate_shortest_distance(data)
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Relating Opening Angle to Distance of Closest Approach',
                    [r'Opening Angle [$\pi$ rad]', 'Closest Distance [m]'])
        ax.grid()
    
    h = ax.hist2d(angles / np.pi, distance, bins = [20,20], weights = data['weight'])
    fig.colorbar(h[3], ax = ax, label = 'Number of Vertices')
    
    if outfile != None:
        plt.savefig(outfile)
        plt.close()
        
    return ax
    
    

def main(filename = None, outfile_prefix = None):
    
    if filename == None:
        filename, outfile_prefix = sys.argv[1:]
        
    mass, mixing = map(float, filename[:-7].split('_')[1:])
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)['Data']
    
    make_opening_histogram(data, filename, outfile_prefix + 'opening_angle.pdf', log = False)
    make_opening_histogram(data, filename, outfile_prefix + 'opeing_angle_log.pdf')
    make_closest_histogram(data, filename, outfile_prefix + 'closest.pdf', log = False)
    make_closest_histogram(data, filename, outfile_prefix + 'closest_log.pdf')
    compare_opening_closest(data,filename, outfile_prefix + 'opening_closest.pdf')
    compare_opening_magnitude(data, filename, outfile_prefix + 'opening_magnitude.pdf')
    
if __name__ == '__main__':
    #main()
    filename = 'sim_0.1_2.335721469090121e-06.pickle'

    with open(filename, 'rb') as f:
        data = pickle.load(f)['Data']
        
    make_opening_histogram(data)
    make_opening_histogram(data, log = False)
    
    make_closest_histogram(data)
    make_closest_histogram(data, log = False)
