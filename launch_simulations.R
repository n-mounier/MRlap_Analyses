#### Check effect of sample overlap on MR estimate ####
################ Using Simulated Data #################

## Note, this script launch 100 simulations using "simulate_data.R" function(s)



#### Main Function - no correlated pleiotropy ####
launch_simulation  <- function(scenario,
                               n_A=20000, 
                               n_B=20000,
                               pi_x=0.001, # proportion of SNPs affecting X
                               h2_x=0.4,  # heritability of X
                               pi_y=0.01, # proportion of SNPs affecting Y
                               h2_y=0.2 , # heritability of Y
                               pi_xy=0, # proportion of X SNPs also affecting Y, uncorrelated pleiotropy
                               kappa_x=0.3, # effect of U on X
                               kappa_y=0.5, # effect of U in Y
                               alpha=0.2){
  library(data.table)
  library(tidyverse)
  options(datatable.fread.datatable=FALSE)
  options("stringsAsFactors" = FALSE)
  library(rslurm)
  
  
  ### JURA or HPC1?
  server = Sys.info()["nodename"]
  if(stringr::str_detect(server, "chuv.vital-it.ch")){
    JURA = TRUE
  } else if(stringr::str_detect(server, ".cluster")){
    JURA = FALSE
  } else {
    stop("problem with server idenfication")
  }
  
  
  if(JURA){
    my_slurm_options = list(partition = "normal", 
                            time="2-00:00:00",  
                            `cpus-per-task`=1, 
                            mem="10G")
  } else {
    my_slurm_options = list(partition = sample(c("cluster", "sgg", "cluster2"), 1), 
                            time="1-00:00:00",  mem="10G",
                            `cpus-per-task`=1)
  }
  
  #### check parameters ####
  # check: maximum number of individuals
  if(n_A + n_B > 370000) stop("Maximum number of individuals (370,000) exceeded.")
  # check: simulations paramaters
  if(h2_x+kappa_x^2>=1) stop("X can not have a variance of 1 if h2_x + kappa_x^2 >= 1.")
  if(h2_y+alpha^2+kappa_y^2+2*alpha*kappa_x*kappa_y>=1) stop("Y can not have a variance of 1 if h2_y + alpha^2 + kappa_y^2 + 2 * alpha * kappa_y * kappa_x >= 1.")
  
  if(pi_xy>pi_x) stop("pi_xy should be smaller than pi_x")
  if(pi_xy>pi_y) stop("pi_xy should be smaller than pi_y")
  
  if(JURA){
    setwd("~/scratch_kuta//nmounier/projects/SampleOverlap/Data/Simulations")
  }  else{
    setwd("/data/sgg3/ninon/projects/SampleOverlap/Data/Simulations/")
  }
  
  # creates data folder
  if(!dir.exists(scenario)) dir.create(scenario)
  setwd(scenario)
  # also creates results folder
  res_folder = paste0("../../../Results/Simulations/", scenario)
  if(!dir.exists(res_folder)) dir.create(res_folder)
  
  # check: should not be any results for these IDs already (use ID_seed.txt)
  #if(any(file.exists(paste0(1:100, "_seed.txt")))) stop("Some IDs already used for this scenario.")
  #source("../../../Scripts/simulate_data.R")
  
  
  # params : to write in folder
  params = data.frame(scenario = scenario,
                      n_A = n_A,
                      n_B = n_B,
                      pi_x = pi_x, # proportion of SNPs affecting X
                      h2_x = h2_x,  # heritability of X
                      pi_y = pi_y, # proportion of SNPs affecting Y
                      h2_y = h2_y,  # heritability of Y
                      pi_xy = pi_xy, # correlated pleiotropy
                      kappa_x = kappa_x, # effect of U on X
                      kappa_y = kappa_y, # effect of U in Y
                      alpha = alpha)
  write.table(params, paste0(res_folder, "/params.csv"), quote=F, row.names=F)
  
  
  # launch simulations 1:100
  id_tolaunch = 1
  launched = 0
  
  while(id_tolaunch<100){
    
    ids = c(id_tolaunch:min(id_tolaunch+5,100))
    id_tolaunch = max(ids) + 1
    
    
    # params for rslurm (add IDs)
    params = list(simulation_IDs = ids,
                  scenario = params$scenario,
                  n_A = params$n_A,
                  n_B = params$n_B,
                  pi_x = params$pi_x,
                  h2_x = params$h2_x,
                  pi_y = params$pi_y,
                  h2_y = params$h2_y,
                  kappa_x = params$kappa_x,
                  kappa_y = params$kappa_y,
                  alpha = params$alpha)
    
    
    slurm_call(simulate_data, params,
               jobname = paste0(scenario, "_", launched),
               # specify where packages are stored
               libPaths=.libPaths(),
               # specify partition (sgg/cluster/cluster2)
               slurm_options = my_slurm_options,
               submit = TRUE)
    
    
    
    launched = launched+1
    
    # launch bunch of 2*6, then wait 2 hours
    if(launched %% 2 == 0)   Sys.sleep(2*60*60)
    #Sys.sleep(4*60*60)
  }
  
  
  # clean
  # should wait until all jobs are done!!!
  all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0(scenario, "_", 1:launched, collapse=",")), intern=TRUE))
  while(length(all_myjobs)>1){
    # test every 5 minutes
    Sys.sleep(5*60)
    all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0(scenario, "_", 1:launched, collapse=",")), intern=TRUE))
  }
  
  # remove all BGENIE outputs
  system("rm *.out.gz")
  
  # remove all _rslurm folders
  system("rm -rf _rslurm*")
  system("rm results_*.RDS")
  
  # that means that in the folders, only these files are left (for each ID):
  # ID_seed.txt                   
  # ID_effectiveSNPs_X.csv 
  # ID_effectiveSNPs_Y.csv
  # GWAS_X_ID.tsv
  # GWAS_Y0_ID.tsv
  # ...
  # GWAS_Y100_ID.tsv
  # + phenofile (one for multiple IDs)
  
} 

#### standard settings ####
# launch_simulation("noY",
#                   n_A=20000, 
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   alpha=0.2)


#### no causal effect ####
# launch_simulation("noY_nocausaleffect",
#                   n_A=20000, 
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   alpha=0)


#### stronger U ####
# launch_simulation("noY_strongerU",
#                   n_A=20000, 
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=0.5, # effect of U on X
#                   kappa_y=0.8, # effect of U in Y
#                   alpha=0.2)


#### weaker U ####
# launch_simulation("noY_weakerU",
#                   n_A=20000, 
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=0.15, # effect of U on X
#                   kappa_y=0.3, # effect of U in Y
#                   alpha=0.2)


#### realistic settings ####
# launch_simulation("noY_realistic",
#                   n_A=100000,
#                   n_B=100000,
#                   pi_x=0.005, # proportion of SNPs affecting X
#                   h2_x=0.2,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   alpha=0.1)


#### negative U ####
# launch_simulation("noY_negativeU",
#                   n_A=20000, 
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=-0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   alpha=0.2)



#### uncorrelated pleiotropy #### 
# launch_simulation("withY_60",
#                   n_A=20000,
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0.002, # proportion of SNPs affecting Y
#                   h2_y=0.3 , # heritability of Y
#                   pi_xy = 60/100 * 0.001,
#                   kappa_x=0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   alpha=0.2)


#### Main Function - correlated pleiotropy ####
launch_simulation_genU  <- function(scenario,
                                    n_A=20000, 
                                    n_B=20000,
                                    pi_x=0.001, # proportion of SNPs affecting X
                                    h2_x=0.4,  # heritability of X
                                    pi_y=0.002, # proportion of SNPs affecting Y
                                    h2_y=0.2 , # heritability of Y
                                    pi_xy=0, # proportion of X SNPs also affecting Y, uncorrelated pleiotropy
                                    kappa_x=0.3, # effect of U_e on X
                                    kappa_y=0.5, # effect of U_e in Y
                                    pi_u=0.005, # proportion of SNPs affecting U_g
                                    h2_u=0.1,  # heritability of U_g
                                    q_x = 0.1, # effect of U_g on X
                                    q_y = 0.2, # effect of U_g on Y
                                    alpha=0.2){
  library(data.table)
  library(tidyverse)
  options(datatable.fread.datatable=FALSE)
  options("stringsAsFactors" = FALSE)
  library(rslurm)
  
  
  ### JURA or HPC1?
  server = Sys.info()["nodename"]
  if(stringr::str_detect(server, "chuv.vital-it.ch")){
    JURA = TRUE
  } else if(stringr::str_detect(server, ".cluster")){
    JURA = FALSE
  } else {
    stop("problem with server idenfication")
  }
  
  
  if(JURA){
    my_slurm_options = list(partition = "normal", 
                            time="2-00:00:00",  
                            `cpus-per-task`=1, 
                            mem="10G")
  } else {
    my_slurm_options = list(partition = sample(c( "sgg", "cluster2"), 1), 
                            time="1-00:00:00",  mem="10G",
                            `cpus-per-task`=1)
  }
  
  #### check parameters ####
  # check: maximum number of individuals
  if(n_A + n_B > 370000) stop("Maximum number of individuals (370,000) exceeded.")
  # check: simulations paramaters
  if(h2_x+kappa_x^2>=1) stop("X can not have a variance of 1 if h2_x + kappa_x^2 >= 1.")
  if(h2_y+alpha^2+kappa_y^2+2*alpha*kappa_x*kappa_y>=1) stop("Y can not have a variance of 1 if h2_y + alpha^2 + kappa_y^2 + 2 * alpha * kappa_y * kappa_x >= 1.")
  
  if(pi_xy>pi_x) stop("pi_xy should be smaller than pi_x")
  if(pi_xy>pi_y) stop("pi_xy should be smaller than pi_y")
  
  if(JURA){
    setwd("~/scratch_kuta//nmounier/projects/SampleOverlap/Data/Simulations")
  }  else{
    setwd("/data/sgg3/ninon/projects/SampleOverlap/Data/Simulations/")
  }
  
  # creates data folder
  if(!dir.exists(scenario)) dir.create(scenario)
  setwd(scenario)
  # also creates results folder
  res_folder = paste0("../../../Results/Simulations/", scenario)
  if(!dir.exists(res_folder)) dir.create(res_folder)
  
  # check: should not be any results for these IDs already (use ID_seed.txt)
  if(any(file.exists(paste0(1:100, "_seed.txt")))) stop("Some IDs already used for this scenario.")
  #source("../../../Scripts/simulate_data.R")
  
  
  # params : to write in folder
  params = data.frame(scenario = scenario,
                      n_A = n_A,
                      n_B = n_B,
                      pi_x = pi_x, # proportion of SNPs affecting X
                      h2_x = h2_x,  # heritability of X
                      pi_y = pi_y, # proportion of SNPs affecting Y
                      h2_y = h2_y,  # heritability of Y
                      pi_xy = pi_xy, # correlated pleiotropy
                      kappa_x = kappa_x, # effect of U on X
                      kappa_y = kappa_y, # effect of U in Y
                      pi_u = pi_u, # proportion of SNPs affecting U_g
                      h2_u= h2_u,  # heritability of U_g
                      q_x = q_x, # effect of U_g on X
                      q_y = q_y, # effect of U_g on Y
                      alpha = alpha)
  write.table(params, paste0(res_folder, "/params.csv"), quote=F, row.names=F)
  
  
  # launch simulations 1:100
  id_tolaunch = 1
  launched = 0
  
  while(id_tolaunch<100){
    
    ids = c(id_tolaunch:min(id_tolaunch+5,100))
    id_tolaunch = max(ids) + 1
    
    
    # params for rslurm (add IDs)
    params = list(simulation_IDs = ids,
                  scenario = params$scenario,
                  n_A = params$n_A,
                  n_B = params$n_B,
                  pi_x = params$pi_x,
                  h2_x = params$h2_x,
                  pi_y = params$pi_y,
                  h2_y = params$h2_y,
                  pi_xy = params$pi_xy, # correlated pleiotropy
                  kappa_x = params$kappa_x,
                  kappa_y = params$kappa_y,
                  pi_u = params$pi_u, # proportion of SNPs affecting U_g
                  h2_u= params$h2_u,  # heritability of U_g
                  q_x = params$q_x, # effect of U_g on X
                  q_y = params$q_y, # effect of U_g on Y
                  alpha = params$alpha)
    
    
    slurm_call(simulate_data_genU, params,
               jobname = paste0(scenario, "_", launched),
               # specify where packages are stored
               libPaths=.libPaths(),
               # specify partition (sgg/cluster/cluster2)
               slurm_options = my_slurm_options,
               submit = TRUE)
    
    
    
    launched = launched+1
    
    # launch bunch of 2*6, then wait 4 hours
    if(launched %% 2 == 0)   Sys.sleep(4*60*60)
    
  }
  
  
  # clean
  # should wait until all jobs are done!!!
  all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0(scenario, "_", 1:launched, collapse=",")), intern=TRUE))
  while(length(all_myjobs)>1){
    # test every 5 minutes
    Sys.sleep(5*60)
    all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0(scenario, "_", 1:launched, collapse=",")), intern=TRUE))
  }
  
  # remove all BGENIE outputs
  system("rm *.out.gz")
  
  # remove all _rslurm folders
  system("rm -rf _rslurm*")
  system("rm results_*.RDS")
  
  # that means that in the folders, only these files are left (for each ID):
  # ID_seed.txt                   
  # ID_effectiveSNPs_X.csv 
  # ID_effectiveSNPs_Y.csv
  # GWAS_X_ID.tsv
  # GWAS_Y0_ID.tsv
  # ...
  # GWAS_Y100_ID.tsv
  # + phenofile (one for multiple IDs)
  
} 


#### strong genetic confounder #### 
# launch_simulation_genU("genU",
#                        n_A=20000,
#                        n_B=20000,
#                        pi_x=0.001, # proportion of SNPs affecting X
#                        h2_x=0.4,  # heritability of X
#                        pi_y=0, # proportion of SNPs affecting Y
#                        h2_y=0 , # heritability of Y
#                        pi_xy = 0, # no uncorrelated pleiotropy
#                        kappa_x=0.3, # effect of U on X
#                        kappa_y=0.5, # effect of U in Y
#                        pi_u=0.0005, # proportion of SNPs affecting U_g
#                        h2_u=0.3,  # heritability of U_g
#                        q_x = 0.5, # effect of U_g on X
#                        q_y = 0.7, # effect of U_g on Y
#                        alpha=0.2)

#### weaker genetic confounder #### 
# launch_simulation_genU("genU_moderate",
#                   n_A=20000,
#                   n_B=20000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   pi_xy = 0, # no uncorrelated pleiotropy
#                   kappa_x=0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   pi_u=0.0001, # proportion of SNPs affecting U_g
#                   h2_u=0.2,  # heritability of U_g
#                   q_x = 0.4, # effect of U_g on X
#                   q_y = 0.3, # effect of U_g on Y
#                   alpha=0.2)



#### Main Function - case control (analysed as continuous) ####
launch_simulation_CC  <- function(scenario,
                                  n_A=20000, 
                                  n_B=20000,
                                  pi_x=0.001, # proportion of SNPs affecting X
                                  h2_x=0.4,  # heritability of X
                                  pi_y=0.01, # proportion of SNPs affecting Y
                                  h2_y=0.2 , # heritability of Y
                                  kappa_x=0.3, # effect of U on X
                                  kappa_y=0.5, # effect of U in Y
                                  alpha=0.2,
                                  prevalence=0.1){
  library(data.table)
  library(tidyverse)
  options(datatable.fread.datatable=FALSE)
  options("stringsAsFactors" = FALSE)
  library(rslurm)
  
  
  ### JURA or HPC1?
  server = Sys.info()["nodename"]
  if(stringr::str_detect(server, "chuv.vital-it.ch")){
    JURA = TRUE
  } else if(stringr::str_detect(server, ".cluster")){
    JURA = FALSE
  } else {
    stop("problem with server idenfication")
  }
  
  
  if(JURA){
    my_slurm_options = list(partition = "normal", 
                            time="2-00:00:00",  
                            `cpus-per-task`=1, 
                            mem="10G")
  } else {
    my_slurm_options = list(partition = sample(c("cluster", "sgg", "cluster2"), 1), 
                            time="1-00:00:00",  mem="10G",
                            `cpus-per-task`=1)
  }
  
  #### check parameters ####
  # check: maximum number of individuals
  if(n_A + n_B > 370000) stop("Maximum number of individuals (370,000) exceeded.")
  # check: simulations paramaters
  if(h2_x+kappa_x^2>=1) stop("X can not have a variance of 1 if h2_x + kappa_x^2 >= 1.")
  if(h2_y+alpha^2+kappa_y^2+2*alpha*kappa_x*kappa_y>=1) stop("Y can not have a variance of 1 if h2_y + alpha^2 + kappa_y^2 + 2 * alpha * kappa_y * kappa_x >= 1.")
  
  
  if(JURA){
    setwd("~/scratch_kuta//nmounier/projects/SampleOverlap/Data/Simulations")
  }  else{
    setwd("/data/sgg3/ninon/projects/SampleOverlap/Data/Simulations/")
  }
  
  # creates data folder
  if(!dir.exists(scenario)) dir.create(scenario)
  setwd(scenario)
  # also creates results folder
  res_folder = paste0("../../../Results/Simulations/", scenario)
  if(!dir.exists(res_folder)) dir.create(res_folder)
  
  # check: should not be any results for these IDs already (use ID_seed.txt)
  if(any(file.exists(paste0(1:100, "_seed.txt")))) stop("Some IDs already used for this scenario.")
  #source("../../../Scripts/simulate_data.R")
  
  
  # params : to write in folder
  params = data.frame(scenario = scenario,
                      n_A = n_A,
                      n_B = n_B,
                      pi_x = pi_x, # proportion of SNPs affecting X
                      h2_x = h2_x,  # heritability of X
                      pi_y = pi_y, # proportion of SNPs affecting Y
                      h2_y = h2_y,  # heritability of Y
                      kappa_x = kappa_x, # effect of U on X
                      kappa_y = kappa_y, # effect of U in Y
                      alpha = alpha,
                      prevalence=prevalence)
  write.table(params, paste0(res_folder, "/params.csv"), quote=F, row.names=F)
  
  
  # launch simulations 1:100
  id_tolaunch = 1
  launched = 0
  
  while(id_tolaunch<100){
    
    ids = c(id_tolaunch:min(id_tolaunch+9,100))
    id_tolaunch = max(ids) + 1
    
    
    # params for rslurm (add IDs)
    params = list(simulation_IDs = ids,
                  scenario = params$scenario,
                  n_A = params$n_A,
                  n_B = params$n_B,
                  pi_x = params$pi_x,
                  h2_x = params$h2_x,
                  pi_y = params$pi_y,
                  h2_y = params$h2_y,
                  kappa_x = params$kappa_x,
                  kappa_y = params$kappa_y,
                  alpha = params$alpha,
                  prevalence = params$prevalence)
    
    
    slurm_call(simulate_data_CC, params,
               jobname = paste0(scenario, "_", launched),
               # specify where packages are stored
               libPaths=.libPaths(),
               # specify partition (sgg/cluster/cluster2)
               slurm_options = my_slurm_options,
               submit = TRUE)
    
    
    
    launched = launched+1
    
    # launch bunch of 2*10, then wait 2.5 hours
    if(launched %% 2 == 0)   Sys.sleep(2*90*60)
    
  }
  
  
  # clean
  # should wait until all jobs are done!!!
  all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0(scenario, "_", 1:launched, collapse=",")), intern=TRUE))
  while(length(all_myjobs)>1){
    # test every 5 minutes
    Sys.sleep(5*60)
    all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0(scenario, "_", 1:launched, collapse=",")), intern=TRUE))
  }
  
  # remove all BGENIE outputs
  system("rm *.out.gz")
  
  # remove all _rslurm folders
  system("rm -rf _rslurm*")
  system("rm results_*.RDS")
  
  # that means that in the folders, only these files are left (for each ID):
  # ID_seed.txt                   
  # ID_effectiveSNPs_X.csv 
  # ID_effectiveSNPs_Y.csv
  # GWAS_X_ID.tsv
  # GWAS_Y0_ID.tsv
  # ...
  # GWAS_Y100_ID.tsv
  # + phenofile (one for multiple IDs)
  
} 


#### case control #### 
# launch_simulation("noY_CC",
#                   n_A=100000,
#                   n_B=100000,
#                   pi_x=0.001, # proportion of SNPs affecting X
#                   h2_x=0.4,  # heritability of X
#                   pi_y=0, # proportion of SNPs affecting Y
#                   h2_y=0 , # heritability of Y
#                   kappa_x=0.3, # effect of U on X
#                   kappa_y=0.5, # effect of U in Y
#                   alpha=0.2,
#                   prevalence=0.1)