# -*- coding: utf-8 -*-
"""
This script makes plots for the desired (flavour, mass, mixing_sq) triples given
"""
import sys
import matplotlib.pyplot as plt

import timeit

import make_plots as mp

sys.path.insert(0,'../../FastSim_Additions/')
from Run_Combiner import main as combine
from Additions import remove_infs

plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6,4), 'legend.frameon':True, 'font.size':14})
    
def heat_maps(data, key):
    mp.compare_opening_closest(data, key, outfile = f'plot_{key[0]}_{key[1]}_{key[2]}_twoangles.pdf')
    mp.compare_opening_magnitude(data, key, outfile = f'plot_{key[0]}_{key[1]}_{key[2]}_magnitude.pdf')
    
def multi_curve(data_dict, labels = [], log = True):
    
    funcs = [mp.make_opening_histogram, mp.make_closest_histogram, mp.make_phi_angle_histogram]
    plot_labels = [[r'Maximum Opening Angle [$\pi$ rad]', 'Number of Vertices', (0,1), (0,0.65)],
                   ['Shortest Distance to CMS IP [m]', 'Number of Vertices', (0,200), (0,0.9)],
                   [r'Reconstruction Angle [$\pi$ rad]', 'Number of Vertices', (0,0.27), (0,0.9)]]
    
    suffix = ['openingangle', 'shortestdist','reconstructionangle']
    
    for plot_func, labels, suff in zip(funcs, plot_labels, suffix):
        fig, ax = plt.subplots(ncols = 2,nrows = 1, figsize = (13,5))
        
        for key, l in zip(data_dict.keys(), ['-','--',':']):
            data = data_dict[key]
            
            plot_func(data, key, ax = ax[0], line = l, norm = True)
            plot_func(data, key, ax = ax[1], line = l)
            
        fig.suptitle(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000$ fb$^{-1}$)')
        for i in [0,1]:
            ax[i].set_xlabel(labels[0])
            ax[1].set_ylabel(labels[1])
            ax[0].set_ylabel('Fraction of Vertices')
            ax[i].grid()
            ax[i].legend(fontsize = 10)
        
        #_, max_y = ax[1].get_ylim()
        #ax[1].set_ylim(1 ,10e8)
        ax[0].set_xlim(*labels[2])
        ax[1].set_xlim(*labels[2])
        ax[0].set_ylim(*labels[3])
        ax[1].set_yscale('log')
        
        outfile = f'plot_poster_{suff}_4.svg'
        plt.savefig(outfile)
        plt.show()
        plt.close()
        

def main(plots_to_make):
    for param_set in plots_to_make:
        data_dict = dict()
        for param in param_set:
            print('Starting Combine')
            start = timeit.default_timer()
            data_w_inf = combine(param)
            new_data = remove_infs(data_w_inf)
            end = timeit.default_timer()
            print('Combine Time: ', end - start)
            #heat_maps(new_data, param)
            data_dict[param] = new_data

        multi_curve(data_dict)
    return data_dict

if __name__ == '__main__':
    plots_to_make = [[('Ue', 0.153993, 1.9306977288832498e-06),
                     ('Ue',1.0838, 3.162277660168379e-07),
                     ('Ue',4.08883,5.1794746792312124e-08)]]
                     #[('Ue', 3.48098000e+00, 1.17876863e-05)]]
    #trial_data = main(plots_to_make)
    main(plots_to_make)
    
    