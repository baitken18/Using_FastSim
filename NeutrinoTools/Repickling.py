# -*- coding: utf-8 -*-
"""
Takes neutrino flux pickles from cern box and repackages into a more favourable form
"""
import sys
import pickle 
import numpy as np

def extract_pickle(filename):
    
    with open(filename, 'rb') as f:
        inputev = pickle.load(f)
        sigmaeff = pickle.load(f)
        weights = pickle.load(f)
        pdgs = pickle.load(f)
        pdgfreqs = pickle.load(f)
        
    return inputev, sigmaeff, weights, pdgs, pdgfreqs

def pdg_histo(pdgs, pdgfreqs):
    return np.array([pdgs,pdgfreqs]).T  

def repackage(data):
    
    return {'inputev': data[0],
            'sigmaeff': data[1],
            'weights': np.array(data[2]),
            'pdgs':pdg_histo(*data[3:])}

def main(infile = None, outfile = None, save = False):
    
    if infile == None:
        infile,outfile = sys.argv[1:]
        
    data = extract_pickle(infile)
    
    new_data = repackage(data)
    
    if save:
        with open(outfile, 'wb') as f:
            pickle.dump(new_data, f)
        
    return new_data
        
if __name__ == '__main__':
    infile = '../../AllSimData/nue_0-1Gev_100kev_gntp.0.gst.pickle'
    outfile = 'test.pickle'
    main(infile,outfile)

