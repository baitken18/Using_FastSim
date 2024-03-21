from Run_Optimizer import get_optimal, make_bgd_hists
import glob

import matplotlib.pyplot as plt
import numpy as np
import pickle 

import sys
sys.path.insert(0, '../../FastSim_Additions/')
from Run_Combiner import main as combine

plt.style.use('science')
plt.rcParams.update({'figure.figsize':(6,4), 'font.size':14, 'legend.frameon':True})

def format_plot(ax):
    ax.grid()
    ax.legend()
    ax.set_title(r'MATHUSLA100 No Walls Geometry ($\int\mathcal{L}dt = 3000 fb^{-1}$)')
    ax.set_xlabel('Signal Requirement (With DV3)')
    ax.set_ylabel('Background Rejection (With DV3)')

def get_scatter_pts(RHN_params):
    RHN_data = combine(RHN_params)
    
    #Bgd_data = make_bgd_hists(glob.glob('../NeutrinoResults/*.pickle'), track_thresh = 3)
    with open('BgdHistDV3.pickle', 'rb') as f: 
        Bgd_data = pickle.load(f)
    print('Data Loaded')
        
    signal_reqs = np.arange(80, 99, 1) / 100 
    bgd_cut = [get_optimal(RHN_data, Bgd_data, con = sig, track_thresh = 3)[1] for sig in signal_reqs]
    print(bgd_cut)
    
    return signal_reqs, bgd_cut

def main():
    #Define list of RHN Params
    all_params = [('Ue',4.08883,5.1794746792312124e-08),('Ue',1.0838, 3.162277660168379e-07)]
    fig, ax = plt.subplots(nrows = 1,ncols = 1)
    
    for pars in all_params:
        x,y = get_scatter_pts(pars)
        ax.scatter(x,y, label = fr'RHN (M={pars[1]}GeV, $\vert U_e\vert^2$ = {pars[2]:.3g})')
    
    format_plot(ax)
    plt.savefig('ROCplotDV3_V2.pdf')
    plt.show()
    plt.close()

if __name__ == '__main__':
    main()

