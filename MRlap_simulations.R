### Launch MRlap & get MRlap results for simulated data ###

if(!utils::packageVersion("MRlap") == "0.0.2.0") stop("for paper, use MRlap v. 0.0.2.0")


launch_MRlap <- function(scenario, my_IDs=1:100, my_thresholds = c(1e-08, 5e-08, 1e-07, 5e-07, 1e-06)){
  
  library(rslurm)
  get_job_status <- function (slr_job) {
    if (!(class(slr_job) == "slurm_job"))
      stop("input must be a slurm_job")
    stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname),
                                    intern = TRUE))
    if (length(stat) > 1) {
      res = "Job running or in queue."
    }
    else {
      res = "Job completed or stopped."
    }
    return(res)
  }
  
  
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
    setwd(paste0("~/scratch_kuta//nmounier/projects/SampleOverlap/Data/Simulations/", scenario))
  }  else{
    setwd(paste0("/data/sgg3/ninon/projects/SampleOverlap/Data/Simulations/", scenario))
  }
  
  
  if(JURA){
    my_slurm_options = list(partition = "normal",
                            time="1-00:00:00",
                            `cpus-per-task`=1,
                            mem="8G")
  } else {
    my_slurm_options = list(partition = "sgg",
                            time="1-00:00:00",
                            `cpus-per-task`=1)
  }
  
  pars = data.frame(my_ID=my_IDs, scenario=scenario)
  
  
  analyse_one <- function(my_ID, scenario){
    library(data.table)
    library(tidyverse)
    options(datatable.fread.datatable=FALSE)
    options("stringsAsFactors" = FALSE)
    library(rslurm)
    library(MRlap)
    
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
      setwd(paste0("~/scratch_kuta//nmounier/projects/SampleOverlap/Data/Simulations/", scenario))
    }  else{
      setwd(paste0("/data/sgg3/ninon/projects/SampleOverlap/Data/Simulations/", scenario))
    }
    
    data_folder = getwd()
    ## format data with sample size
    X = data.table::fread(paste0("GWAS_X_", my_ID, ".tsv"))
    res_folder = paste0("../../../Results/Simulations/", scenario)
    params = data.table::fread(paste0(res_folder, "/params.csv"))
    X$N = params$n_A
    
    setwd(res_folder)
    dir.create(paste0("MRlap_", my_ID))
    setwd(paste0("MRlap_", my_ID))

    
    ## ld data
    if(!file.exists("eur_w_ld_chr")){
      if(JURA){
        system("ln -s ~/scratch_kuta/nmounier/data/LDscores/eur_w_ld_chr/")
      } else {
        system("ln -s /data/sgg2/aaron/shared/ldsc_nice_script/ld_scores_precomputed/eur_w_ld_chr/")
      }
    }
    
    if(!file.exists("w_hm3.noMHC.snplist")){
      if(JURA){
        system("ln -s ~/scratch_kuta/nmounier/data/GenomicSEM/w_hm3.noMHC.snplist")
      } else {
        system("ln -s /data/sgg2/ninon/projects/Software_Tests/Test_GenomicSEM//Lifespan/Data/w_hm3.noMHC.snplist")
      }
    }
    
    for(my_overlap in c(0,25,50,75,100)){
      Y = data.table::fread(paste0(data_folder, "/GWAS_Y", my_overlap, "_", my_ID, ".tsv"))
      Y$N =  params$n_B
      
      for(my_threshold in my_thresholds){
        res = MRlap(exposure=X,
                    exposure_name = "exp",
                    outcome=Y,
                    outcome_name = "out",
                    ld="eur_w_ld_chr",
                    hm3="w_hm3.noMHC.snplist",
                    MR_threshold = my_threshold,
                    MR_pruning_dist = 500,
                    MR_pruning_LD = 0,
                    MR_reverse = 1e-3,
                    save_logfiles = FALSE,
                    verbose = TRUE)
        # save RDS file
        saveRDS(res, paste0("res_", my_overlap, "_", my_threshold, ".RDS"))
      }
      
    }
    
    
    
  }
  
  
  sjob_scenario <- slurm_apply(analyse_one, pars,
                               jobname = paste0('analysis_', scenario),
                               nodes = nrow(pars), cpus_per_node = 1,
                               # specify where packages are stored
                               #libPaths=.libPaths(),
                               # specify partition (sgg/cluster/cluster2)
                               slurm_options = my_slurm_options,
                               submit = TRUE)
  
  wait=T
  while(wait){
    # try once every minute, that's enough
    Sys.sleep(60)
    #print(get_job_status(sjob_getEffects))
    if(get_job_status(sjob_scenario) == "Job completed or stopped.") wait = F
  }
  
}

# launch_MRlap("noY")
# launch_MRlap("noY_negativeU")
# launch_MRlap("noY_nocausaleffect")
# launch_MRlap("noY_strongerU")
# launch_MRlap("noY_weakerU")
# launch_MRlap("noY_realistic")
# launch_MRlap("withY_60")

# launch_MRlap("genU")
# launch_MRlap("genU_moderate")

# launch_MRlap("noY_CC")



get_MRlap_results <- function(scenario, my_IDs=1:100, my_thresholds=c(1e-8, 5e-8, 1e-7, 5e-7, 1e-6), suffix=""){
  
  server = Sys.info()["nodename"]
  if(stringr::str_detect(server, "chuv.vital-it.ch")){
    JURA = TRUE
  } else if(stringr::str_detect(server, ".cluster")){
    JURA = FALSE
  } else {
    stop("problem with server idenfication")
  }
  
  
  if(JURA){
    setwd(paste0("~/scratch_kuta//nmounier/projects/SampleOverlap/Data/Simulations/", scenario))
  }  else{
    setwd(paste0("/data/sgg3/ninon/projects/SampleOverlap/Data/Simulations/", scenario))
  }
  
  
  res_folder = paste0("../../../Results/Simulations/", scenario)
  setwd(res_folder)
  
  res = expand.grid(ID=my_IDs, overlap=c(0,25,50,75,100), threshold=my_thresholds)
  
  res$obs_alpha = NA_real_
  res$obs_alpha_se = NA_real_
  res$corrected_alpha = NA_real_
  res$corrected_alpha_se = NA_real_
  res$m = NA_real_
  
  for(my_ID in my_IDs){
    print(my_ID)
    setwd(paste0("MRlap_", my_ID))
    for(my_overlap in c(0,25,50,75,100)){
      for(my_threshold in my_thresholds){
        MRlap_res = readRDS(paste0("res_", my_overlap, "_", my_threshold, ".RDS"))
        res[res$ID == my_ID & res$overlap == my_overlap & res$threshold == my_threshold, "m"] = MRlap_res$MRcorrection$m_IVs
        res[res$ID == my_ID & res$overlap == my_overlap & res$threshold == my_threshold, "obs_alpha"] = MRlap_res$MRcorrection$observed_effect
        res[res$ID == my_ID & res$overlap == my_overlap & res$threshold == my_threshold, "obs_alpha_se"] = MRlap_res$MRcorrection$observed_effect_se
        res[res$ID == my_ID & res$overlap == my_overlap & res$threshold == my_threshold, "corrected_alpha"] = MRlap_res$MRcorrection$corrected_effect
        res[res$ID == my_ID & res$overlap == my_overlap & res$threshold == my_threshold, "corrected_alpha_se"] = MRlap_res$MRcorrection$corrected_effect_se
        
      }
    }
    setwd("..")
  }
  write.table(res, paste0("MRlap_corrected", suffix, ".csv"), sep=",", quote=F, row.names=F)
  
  
}

# get_MRlap_results("noY")
# get_MRlap_results("noY_negativeU")
# get_MRlap_results("noY_nocausaleffect")
# get_MRlap_results("noY_strongerU")
# get_MRlap_results("noY_weakerU")
# get_MRlap_results("noY_realistic")
# get_MRlap_results("withY_60")

# get_MRlap_results("genU")
# get_MRlap_results("genU_moderate")

# get_MRlap_results("noY_CC")
