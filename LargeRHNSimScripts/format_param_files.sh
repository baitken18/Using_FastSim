for file in ~/Documents/5thYear/Phys499/AllSimData/MATHUSLA_LLPfile_RHN_"#1"/RHN_"#1"_hadronic_decays_geant/*
do
python ~/Documents/5thYear/Phys499/FastSim_Additions/Geant_to_hadron_style.py $file
done
python ~/Documents/5thYear/Phys499/Using_FastSim/LargeRHNSimScripts/param_file_writer.py "#1"
python ~/Documents/5thYear/Phys499/Using_FastSim/LargeRHNSimScripts/multi_run_simulation.py "#2" param_file_"#1".txt
#Potentailly make new directory that dumps the data into it
