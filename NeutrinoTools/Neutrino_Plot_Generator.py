# -*- coding: utf-8 -*-
"""
This script formats all of the plots of the neutrino data
"""

import pickle
import Neutrino_Plots as neu 
import glob 

import matplotlib.pyplot as plt

plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6,4), 'legend.frameon':True, 'font.size':14,
                     'legend.fontsize':12, 'legend.handlelength': 1})

def format_2d(ax):
    ax.set_xlabel(r'Opening Angle [$\pi$ rad]')
    ax.set_ylabel('Shortest distance')
    ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000 fb^{-1}$)')
    
def format_1d(ax, style, log = False):
    ax.set_xlabel(style)
    ax.set_ylabel('Number of Vertices')
    ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000 fb^{-1}$)')
    ax.grid()
    if log: ax.set_yscale('log')
    
def plot_opening(data_sets):
    
    fig, ax = plt.subplots(1,1)
    
    for params in data_sets.keys():
        
        neu.add_opening_curve(data_sets[params], params, ax)
    
    format_1d(ax, r'Opening Angle [$\pi$ rad]')
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.13)
    ax.legend()
    
    plt.savefig(f'../NeutrinoResults/Opening_Backgrounds_{params[0]}.pdf')
    plt.show()
    plt.close()
    
def plot_recon(data_sets):
    
    fig, ax = plt.subplots(1,1)
    
    for params in data_sets.keys():
        
        neu.add_reconstruction_curve(data_sets[params], params, ax)
    
    format_1d(ax, r'Reconstruction Angle [$\pi$ rad]')
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.07)
    ax.legend()
    
    plt.savefig(f'../NeutrinoResults/Recon_Backgrounds_{params[0]}.pdf')
    plt.show()
    plt.close()
    
def plot_closest(data_sets):
    
    fig, ax = plt.subplots(1,1)
    
    for params in data_sets.keys():
        
        neu.add_closest_curve(data_sets[params], params, ax)
    
    format_1d(ax, r'Shortest Distance to IP [m]')
    ax.set_xlim(0,200)
    ax.set_ylim(0,0.07)
    ax.legend()
    
    plt.savefig(f'../NeutrinoResults/Closest_Backgrounds_{params[0]}.pdf')
    plt.show()    
    plt.close()
    
def plot_cumulative(filenames):
    
    data_sets = {}
    for file in filenames:
        with open(file, 'rb') as f:
            data_sets[extract_params(file)] = pickle.load(f)
    
    fig, ax = plt.subplots(1,1)
    neu.cumulative_opening_curve(data_sets, ax)
    format_1d(ax, r'Opening Angle [$\pi$ rad]')
    ax.set_xlim(0,1)
    plt.savefig('../NeutrinoResults/Opening_Backgrounds_cumulative.pdf')
    plt.show()    
    plt.close()
    
    fig, ax = plt.subplots(1,1)
    neu.cumulative_closest_curve(data_sets, ax)
    format_1d(ax, r'Shortest Distance to IP [m]')
    plt.savefig('../NeutrinoResults/Closest_Backgrounds_cumulative.pdf')
    plt.show()    
    plt.close()
    
    fig, ax = plt.subplots(1,1)
    neu.cumulative_reconstruction_curve(data_sets, ax)
    ax.set_ylim(0,0.26)
    format_1d(ax, r'Reconstruction Angle [$\pi$ rad]')
    plt.savefig('../NeutrinoResults/Recon_Backgrounds_cumulative.pdf')
    plt.show()    
    plt.close()
    
def extract_params(filename):
    parts = filename.split('\\')[-1].split('_')
    return (parts[0], parts[1][0], parts[1][2:])
    
def main(filenames):
    data_sets = {}
    for file in filenames:
        with open(file, 'rb') as f:
            data_sets[extract_params(file)] = pickle.load(f)
    
    plot_opening(data_sets)
    plot_recon(data_sets)
    plot_closest(data_sets)
    
if __name__ == '__main__':
    '''
    filenames = glob.glob('../NeutrinoResults/n*.pickle')
    plot_cumulative(filenames)
    '''
    
    #filenames = glob.glob('../NeutrinoResults/nue_*.pickle')
    #main(filenames)
    
    #filenames = glob.glob('../NeutrinoResults/nuebar_*.pickle')
    #main(filenames)
    
    filenames = glob.glob('../NeutrinoResults/numu_*.pickle')
    main(filenames)
    
    #filenames = glob.glob('../NeutrinoResults/numubar_*.pickle')
    #main(filenames)
    
            
    
    
    