for file in ~/Documents/5thYear/Phys499/AllSimData/MATHUSLA_LLPfiles_RHN_$1/RHN_$1_hadronic_decays_geant/*
do
python ~/Documents/5thYear/Phys499/FastSim_Additions/Geant_to_hadron_style.py $file
done
echo Geant_formatted populated
python ~/Documents/5thYear/Phys499/Using_FastSim/LargeRHNSimScripts/param_file_writer.py $1
echo Param_file Written
python ~/Documents/5thYear/Phys499/Using_FastSim/LargeRHNSimScripts/multi_process_simulation.py $2 param_file_$1.txt
mv sim*.pickle Finished_Sim_$1/
echo Sim Data Moved 
