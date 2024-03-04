# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:43:14 2024

@author: baitk
"""

import numpy as np
import awkward as ak
import uproot 

def pull_from_root(filename, tree = 'gst', library = 'np'):
    file = uproot.open(f'{filename}:{tree}')
    data = file.arrays(library = library)
    
    return data

def transpose_daughters(E, px, py, pz, pid):
    temp = ak.Array([E, px, py, pz, pid])
    
    builder = ak.ArrayBuilder()
    for i in range(len(E)):
        builder.append(temp[:,i].to_numpy().T)
        
    return builder.snapshot()

def append_new_event(new, saved):
    saved['Ev'] = np.append(saved['Ev'], np.array(new['Ev']))
    saved['weight'] = np.append(saved['weight'], np.array(new['weight']))
    saved['position'] = np.append(saved['position'], new['position'])
    saved['daughters'].append(new['daughters'])
    
    return saved

def initiate_structure():
    
    return {'Ev':np.array([]),
            'weight':np.array([]),
            'position':np.array([]),
            'daughters':ak.ArrayBuilder()}

def freeze_structure(saved):
    saved['position'].reshape(len(saved['weight']), 3)
    saved['daughters'] = saved['daughters'].snapshot()
    
    return saved