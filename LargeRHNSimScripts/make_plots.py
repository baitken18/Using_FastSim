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

import matplotlib.pyplot as plt
plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6,4), 'legend.frameon':True, 'font.size':14})


sys.path.insert(0,'../MomentumExtrapolation/')
from Opening_Angle_Histogram import calculate_opening_angle
from Closest_Approach_Histogram import calculate_shortest_distance

def format_plot(ax, title, axis_labels):
    ax.set_title(title)
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.grid()

def make_opening_histogram(data, outfile = None, ax = None):
    angles, weight_bool = calculate_opening_angle(data)
    weights = data['weight'][weight_bool]
    
    #Numpy computes the histogram.  A bin of 0 is added to both ends so the plot doesn't look funny
    h = list(np.histogram(angles, bins = np.arange(0, np.pi, 0.1), weights = weights))
    h[0] = np.append(np.append(np.zeros(1), h[0]), np.zeros(1))
    h[1] = np.append(np.zeros(1), h[1])
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Opening Angle', [r'Maximum Opening Angle [$\pi$ rad]', 'Number of Vertices'])
        ax.set_yscale('log')
    
    #Change to shift errorbars instead
    l = ax.step(h[1]/np.pi, h[0], where = 'post')
    ax.errorbar((h[1][1:-1] + np.diff(h[1])[1]/2)/np.pi, h[0][1:-1], yerr = np.sqrt(h[0][1:-1]), fmt = '.', color = l[0].get_color())
    
    if outfile != None:
        plt.savefig(outfile)
        plt.close()
    
    return ax

def make_closest_histogram(data, outfile = None, ax = None):
    distances, reconstructable_bool = calculate_shortest_distance(data)
    weights = data['weight'][reconstructable_bool.to_numpy()]
    
    h = list(np.histogram(distances, bins = np.arange(0, np.max(distances), 5), weights = weights))
    h[0] = np.append(np.append(np.zeros(1), h[0]), np.zeros(1))
    h[1] = np.append(np.zeros(1), h[1])
    
    if ax == None:
        fig, ax = plt.subplots(1,1)
        format_plot(ax, 'Distance of Closest Approach', ['Shortest Distance [m]', 'Number of Vertices'])
        ax.set_yscale('log')
    
    l = ax.step(h[1], h[0], where = 'post')
    ax.errorbar(h[1][1:-1] + np.diff(h[1])[2]/2, h[0][1:-1], yerr = np.sqrt(h[0][1:-1]), fmt = '.', color = l[0].get_color())
    
    if outfile != None:
        plt.savefig(outfile)
        
        plt.close()
        
    return ax

def make_phi_angle_histogram(data, outfile = None, ax = None):
    pass

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
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)['Data']
    
    make_opening_histogram(data, outfile_prefix + 'opening_angle.pdf')
    make_closest_histogram(data, outfile_prefix + 'closest.pdf')
    compare_opening_closest(data, outfile_prefix + 'opening_closest.pdf')
    compare_opening_magnitude(data, outfile_prefix + 'opening_magnitude.pdf')
    
if __name__ == '__main__':
    main()

