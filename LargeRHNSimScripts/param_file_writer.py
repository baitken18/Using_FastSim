# -*- coding: utf-8 -*-
"""
This script makes the csv param file
"""
import sys
from glob import glob

def get_length_file(flavour):

    return f'../../AllSimData/MATHUSLA_LLPfiles_RHN_{flavour}/RHNctau{flavour}.dat'

def extract_mass(filename):
    return filename.split('_')[-1][:-4]

def fv_and_prod(flavour):
    
    pairs = []
    
    fv_dir = f'../../AllSimData/MATHUSLA_LLPfiles_RHN_{flavour}/All_RHN_{flavour}_LLPweight4vectors_from_BDWZtau/*'

    fv_files = glob(fv_dir)

    
    prod_dir = f'../../AllSimData/MATHUSLA_LLPfiles_RHN_{flavour}/RHN_{flavour}_hadronic_decays_geant_formatted/*'
    prod_files = glob(prod_dir)
    
    for fv in fv_files:
        for prod in prod_files:
            
            if extract_mass(fv) == extract_mass(prod):
                pairs.append([fv,prod])
    
    return pairs

def main(flavour = None):
    #Change this to use flavour instead
    if flavour == None:
        _, flavour = sys.argv
        
    length_file = get_length_file(flavour)
    
    fv_prod_pairs = fv_and_prod(flavour)
    
    outfile = f'param_file_{flavour}.txt'
    with open(outfile, 'w') as f:
    
        for pair in fv_prod_pairs:
            line = f'{pair[0]} {pair[1]} {length_file} (-12,-1) -1\n'
            f.write(line)

if __name__ == '__main__':
    main()

