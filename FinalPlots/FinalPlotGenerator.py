import sys
import pickle
import glob
from timeit import default_timer

import numpy as np
import awkward as ak

import matplotlib.pyplot as plt

plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6 * 1.5,4 * 1.5), 'legend.frameon':True, 
                     'legend.fontsize':8, 'legend.handlelength': 1})

sys.path.insert(0,'../../FastSim_Additions/')
from Run_Combiner import main as combine

sys.path.insert(0, '../NeutrinoTools/')
import Neutrino_Plots as neu 

def format_1d(ax, xlab, log = False):
    ax.set_xlabel(xlab)
    ax.set_ylabel('Fraction of Vertices')
    ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000 fb^{-1}$)')
    ax.grid()
    if log: ax.set_yscale('log')

def collect_RHN_runs(run_params):
    
    storage = dict([(params, combine(params)) for params in run_params])
    return storage 

def collect_Bgd_runs():
    
    files = glob.glob('../NeutrinoResults/*.pickle')
    storage = dict()
    for i,file in enumerate(files):
        with open(file, 'rb') as f:
            data = pickle.load(f)
        
        storage[i] = data
        
    return storage
    

def make_opening_plot(RHN_runs, Bgd_Runs):
    fig, ax = plt.subplots(1,1)
    format_1d(ax, r'Opening Angle [$\pi$ rad]')
    
    for key in RHN_runs.keys():
        neu.add_opening_curve(RHN_runs[key], key, ax)
        
    neu.cumulative_opening_curve(Bgd_Runs, ax)
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.25)
    ax.legend()
    
    #plt.savefig('Final_Opening_Plot.pdf')
    plt.show()
    plt.close()
    
def make_reconstruction_plot(RHN_runs, Bgd_Runs):
    fig, ax = plt.subplots(1,1)
    format_1d(ax, r'Reconstruction Angle [$\pi$ rad]')
    
    for key in RHN_runs.keys():
        neu.add_reconstruction_curve(RHN_runs[key], key, ax)
        
    neu.cumulative_reconstruction_curve(Bgd_Runs, ax)
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.3)
    
    ax.legend()
    
    #plt.savefig('Final_Recon_Plot.pdf')
    plt.show()
    plt.close()
    
def make_shortest_plot(RHN_runs, Bgd_Runs):
    fig, ax = plt.subplots(1,1)
    format_1d(ax, r'Closest Distance to CMS IP [m]')
    
    for key in RHN_runs.keys():
        neu.add_closest_curve(RHN_runs[key], key, ax)
        
    neu.cumulative_closest_curve(Bgd_Runs, ax)
    
    ax.legend()
    ax.set_xlim(0,200)
    ax.set_ylim(0,0.3)
    
    #plt.savefig('Final_shortest_Plot.pdf')
    plt.show()
    plt.close()
        
    

def main():
    start = default_timer()
    RHN_params = [('Ue',4.08883,5.1794746792312124e-08)] #('Ue',1.0838, 3.162277660168379e-07),
    RHN_data = collect_RHN_runs(RHN_params)
    print('Right-Handed Data Formatted')
    end = default_timer()
    print(f'{end - start} s')
    
    start = default_timer()
    Bgd_data = collect_Bgd_runs()
    print('Background Data Formatted')
    end = default_timer()
    print(f'{end - start} s')
    
    start = default_timer()
    make_opening_plot(RHN_data, Bgd_data)
    print('Finished Opening')
    end = default_timer()
    print(f'{end - start} s')
    
    start = default_timer()
    make_reconstruction_plot(RHN_data, Bgd_data)
    print('Finished Recon')
    end = default_timer()
    print(f'{end - start} s')
    
    start = default_timer()
    make_shortest_plot(RHN_data, Bgd_data)
    print('Finished Short')
    end = default_timer()
    print(f'{end - start} s')
    
    return RHN_data, Bgd_data

if __name__ == '__main__':
    all_data = main()