# -*- coding: utf-8 -*-
"""
This script makes the plot denoting point in Phase space that preserve 
a certain percentage of the signal 
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as mplm
import pickle

from CutAlgorithm import main as get_cuts

sys.path.insert(0,'../../FastSim_Additions')
from Run_Combiner import main as combine 

plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6.5 * 3,4.5), 'legend.frameon':True, 'font.size':14,
                     'legend.fontsize':12, 'legend.handlelength': 1})

def format_entire_plot(fig, ax):
    for axis in ax:
        axis.grid()
        #axis.legend()
        
    fig.suptitle(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000 fb^{-1}$) with DV3 Triggering')
    ax[0].set_xlabel(r'Opening Angle Cut [$\pi$ rad]')
    ax[1].set_xlabel(r'Opening Angle Cut [$\pi$ rad]')
    ax[0].set_ylabel(r'Reconstruction Angle Cut [$\pi$ rad]')
    ax[2].set_xlabel(r'Reconstruction Angle Cut [$\pi$ rad]')
    ax[1].set_ylabel(r'Shortest Distance Cut [m]')
    ax[2].set_ylabel(r'Shortest Distance Cut [m]')
    
    ax[0].set_xlim(0,1)
    ax[0].set_ylim(0,1)
    ax[1].set_xlim(0,1)
    ax[1].set_ylim(0,200)
    ax[2].set_xlim(0,1)
    ax[2].set_ylim(0,200)

def main():
    
    RHN_params = [('Ue',4.08883,5.1794746792312124e-08)]
    
    fig, ax = plt.subplots(nrows = 1,ncols = 3)
    for run in RHN_params:
        data = combine(run)
        opening, recon, shortest = get_cuts(data, track_thresh=3).T
        
        bin_vals = [np.arange(0,1.1,0.1),np.arange(0,1.1, 0.1), np.arange(0,201, 5)]
        cmap = mplm.plasma
        norm = mplc.BoundaryNorm([0,0.5,1], cmap.N)
        
        h0 = ax[0].hist2d(opening/np.pi, recon/np.pi, bins = [bin_vals[0],bin_vals[1]], norm = norm, cmap = cmap)
        h1 = ax[1].hist2d(opening/np.pi, shortest, bins = [bin_vals[0],bin_vals[2]],  norm = norm, cmap = cmap)
        h2 = ax[2].hist2d(recon/np.pi, shortest, bins = [bin_vals[1], bin_vals[2]], norm = norm, cmap = cmap)
        cbar = fig.colorbar(h0[3], cmap = cmap, ticks = [0.25,0.75])
        cbar.ax.set_yticklabels(['Invalid', 'Valid'])
        cbar = fig.colorbar(h1[3], cmap = cmap, ticks = [0.25,0.75])
        cbar.ax.set_yticklabels(['Invalid', 'Valid'])
        cbar = fig.colorbar(h2[3], cmap = cmap, ticks = [0.25,0.75])
        cbar.ax.set_yticklabels(['Invalid', 'Valid'])
        
        
        
    format_entire_plot(fig, ax)
    plt.savefig('CutPhaseSpaceDV3_nocap_2dhist.pdf')
   # with open('CutPhaseSpaceDV2_no_cap_info.pickle', 'wb') as f:
    #    pickle.dump([fig,ax],f)
    
    plt.show()
    plt.close()

if __name__ == '__main__':
    main()