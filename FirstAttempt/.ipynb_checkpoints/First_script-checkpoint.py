# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:28:31 2023

@author: baitk
"""

import os
import sys

import numpy as np

import DetectorSimulation.Detector as Detector
import DetectorSimulation.Vertex as Vertex
import DetectorSimulation.Particle as Particle
import DetectorSimulation.lorentz_transformation as lt
import DetectorSimulation.llp_gun as llp_gun
import DetectorSimulation.llp_gun_new as lg
from Helpers import *

#Initiate the Detector
def initiate_detector():
    param_path = os.path.join(os.getcwd(), '../MATHULSA_FastSim/param_card_CDR.txt')
    detector_benchmark = Detector.Detector(param_path)
    detector_benchmark.clear_detector()
    angles = get_detector_angles(detector_benchmark)
    return detector_benchmark, angles

def inital_data(filename, num_particles):
    LLP_file_path = os.path.join(os.getcwd(), '../Data/{filename}')
    vecs = read_vectors(LLP_file_path, num_particles)
    return vecs

def do_decay(LLP, lifetime, detector_benchmark, angles):
    #Angles defining the detector
    phi_min, phi_max, theta_min, theta_max = angles
    
    LLP_theta = get_theta(LLP[2:])
    
    if (LLP_theta < theta_max) and (LLP_theta > theta_min):
        
        fv_rotated = deal_with_phi(LLP[1:])
        
        pack = get_weight(fv_rotated, LLP[0], lifetime, detector_benchmark)
        
        detector_benchmark.clear_detector()
        
        vert = lg.get_llp('Leptonic2body', LLP[0], pack[1], pack[2], [13,13])
        
        if vert is not None:
            
            detector_benchmark.new_vertex_event(vert)
            
            return detector_benchmark
        
    return None

def simulation_output(detector_benchmark, outfile):
    outfile = join(os.getcwd(), outfile)
    detector_benchmark.detector_display(outfile, show = True)
                
    print('DV2 Criteria:', bool(detector_benchmark.vertex_reconstructed(recon_criteria = 'DVmedium2')))
    print('DV3 Criteria:', bool(detector_benchmark.vertex_reconstructed(recon_criteria = 'DVmedium3')))
    print('LLP trigger:', bool(detector_benchmark.event_pass_trigger()))

def run_simulation(vectors, lifetime, detector_benchmark):
    
    for vec in vectors:
        detector_benchmark = do_decay(vec, lifetime, detector_benchmark)
        simulation_output(detector_benchmark)
    

def main(infile, outfile, num_particles):
    detector_benchmark = initiate_detector()
    vectors = initial_data(data_file, num_particles)
    run_simulation(vectors, lifetime, detector_benchmark)
if __name__ == '__main__':
    try:
        _, infile = sys.argv
        main(infile)
    except ValueError:
        print('Manually Call main()')



