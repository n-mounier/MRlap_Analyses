
# MRlap - Analyses

Here we provide some of the scripts used for the analyses in our paper, "Bias correction for inverse variance weighting Mendelian randomization" [REFERENCE]. It is important to note that these scripts are relying on particular data and software (UK Biobank, BEGENIE) and that they have been designed to be used on specific high performance clusters (folder structure, use of the slurm workload manager, etc). Therefore, they would need to be updated to be used in different settings.   


** simulate_data.R ** : script to simulate the data (required to use launch_simulations.R)   

** launch_simulations.R ** : main script to launch the simulations    

** MRlap_simulations.R ** : main script to analyse simulated data     

** sampling_BMI-BMI.R ** : main script to launch the samplings for BMI-BMI analyses (can be modified to analyse other trait pairs)   

** MRlap_realdata.R ** : main script to analyse real data samplings    
