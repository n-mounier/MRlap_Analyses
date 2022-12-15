#### Check effect of sample overlap on MR estimate ####
### BMI - BMI

## Note, this script launch 100 random samplings using "simulate_data() function


#### Simulations ####
simulate_data <- function(simulation_IDs){
  library(data.table)
  library(tidyverse)
  options(datatable.fread.datatable=FALSE)
  options("stringsAsFactors" = FALSE)
  library(rslurm)
  
  # Use UKBB Data
  # Traits : BMI / BMI
  
  
  # For n repetition
  # get a 100K sample and run BMI GWAS
  # get 9 100K sample and run BMI GWAS
  # with respectively, 100, 95, 90, 75, 50, 25, 10, 5, 0 % of sample overlap with first BMI sample

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
    setwd("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/BMI-BMI2")
  } else {
    setwd("/data/sgg3/ninon/projects/SampleOverlap/Data/BMI-BMI2")
  }
  
  # define functions
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
  
  launch_bgenie <- function(chr, phenofile, folder, name="", JURA=FALSE){
    cat(folder)
    setwd(folder)
    
    if(JURA){
      system(paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/bin/bgenie_v1.3/bgenie_v1.3_static1 ",
                    "--bgen /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", chr, ".bgen ",
                    "--pheno ", phenofile, " ",
                    "--thread 8 ",
                    "--pvals --out ", name, "chr", chr, ".out"))
    } else {
      system(paste0("/data/sgg3/jonathan/bgenie_v1.3/bgenie_v1.3_static1 ",
                    "--bgen /data/sgg3/ninon/projects/SampleOverlap/Data/bgen_subsets/chr", chr, ".bgen ",
                    "--pheno ", phenofile, " ",
                    "--thread 8 ",
                    "--pvals --out ", name, "chr", chr, ".out"))
    }
    
    # Eleonora's command line
    # /data/sgg3/jonathan/bgenie_v1.3/bgenie_v1.3_static1 --bgen
    # /data/sgg3/eleonora/projects/UKBB_GWAS/UK10K_SNPrs/CHR11/chr11.bgen
    # --pheno ../phenofile --pvals --out chr11.out
  }
  
  
  
  # Here, do it differently,
  # directly add 6 (simu) * N (overlap) in a single phenofile 
  
  
  getPhenoFile_OneRepetition <- function(ID){
    seed = round(runif(1, 1, 10^5))
    print(ID)
    print(paste0("seed: ", seed, "\n"))
    set.seed(seed)
    write.table(seed, paste0(ID, "_seed.txt"), sep=",", quote=F, row.names=F, col.names = F)
    
    
    # 1) get samples
    if(JURA){
      bgen_order = data.table::fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
      Data = data.table::fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/SelectedIndividuals_Covariates")
    } else {
      bgen_order = data.table::fread("/data/sgg3/data/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
      Data = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Data/Covariates/SelectedIndividuals_Covariates")
    }
    colnames(Data) = c("ID", "Sex", "Age", paste0("PC", seq(1:40)))
    
    
    
    IDs = Data$ID
    
    IDs_BMI = sample(IDs, 100000)
    IDs_noBMI = IDs[!IDs %in% IDs_BMI]
    # 100, 95, 75, 50, 25, 5, 0
    IDs_100 = IDs_BMI
    IDs_95 = c(sample(IDs_BMI, 95000), sample(IDs_noBMI, 5000))
    IDs_75 = c(sample(IDs_BMI, 75000), sample(IDs_noBMI, 25000))
    IDs_50 = c(sample(IDs_BMI, 50000), sample(IDs_noBMI, 50000))
    IDs_25 = c(sample(IDs_BMI, 25000), sample(IDs_noBMI, 75000))
    IDs_5 = c(sample(IDs_BMI, 5000), sample(IDs_noBMI, 95000))
    IDs_0 = sample(IDs_noBMI, 100000)
    
    # 2) run all GWASs
    if(JURA){     
      Pheno_BMI = data.table::fread(file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno//ukb21067.csv", 
                                    select = c(1, 2202))
    } else {
      Pheno_BMI = data.table::fread(file = "/data/sgg3/data/UKBB/pheno//ukb21067.csv", 
                                    select = c(1, 2202))
    }
    
    Pheno_BMI %>%
      slice(match(IDs_BMI, eid)) -> my_pheno_BMI
    Data %>%
      slice(match(my_pheno_BMI$eid, ID)) -> my_data_BMI
    
    my_data_BMI %>%
      mutate(BMI = my_pheno_BMI[,2]) -> Data_GWAS_BMI
    
    # 1) normalise
    Data_GWAS_BMI %>%
      mutate(normBMI =  qnorm((rank(BMI,na.last="keep")-0.5)/sum(!is.na(BMI)))) -> Data_GWAS_BMI
    
    # 2) adjust
    Pheno_res_BMI <- residuals(lm(paste0("normBMI~Sex+Age+I(Age^2)+", paste0("PC", 1:40, collapse="+")),
                                  data=Data_GWAS_BMI, na.action=na.exclude))
    
    ord = match(bgen_order$ID_1[-1], my_data_BMI$ID)
    
    
    # create pheno file
    bgen_order %>%
      slice(-1) %>%
      transmute(ID=ID_1,
                BMI = replace_na(Pheno_res_BMI[ord], -999)) %>%
      mutate(ID=NULL) -> Phenofile
    
    
    
    # BMI 2
    if(JURA){
      Pheno_BMI2 = data.table::fread(file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb21067.csv", 
                                    select = c(1, 2202))
    } else {
      Pheno_BMI2 = data.table::fread(file = "/data/sgg3/data/UKBB/pheno//ukb21067.csv", 
                                    select = c(1, 2202))
    }
    for(overlap in c(100,95,75,50,25,5,0)){
      
      my_IDs = get(paste0("IDs_", overlap))
      Pheno_BMI2 %>%
        slice(match(my_IDs, eid)) -> my_pheno_BMI2
      Data %>%
        slice(match(my_pheno_BMI2$eid, ID)) -> my_data_BMI2
      
      my_data_BMI2 %>%
        mutate(BMI2 = my_pheno_BMI2[,2]) -> Data_GWAS_BMI2
      
      # 1) normalise
      Data_GWAS_BMI2 %>%
        mutate(normBMI2 =  qnorm((rank(BMI2,na.last="keep")-0.5)/sum(!is.na(BMI2)))) -> Data_GWAS_BMI2
      
      # 2) adjust
      Pheno_res_BMI2 <- residuals(lm(paste0("normBMI2~Sex+Age+I(Age^2)+", paste0("PC", 1:40, collapse="+")),
                                    data=Data_GWAS_BMI2, na.action=na.exclude))
      
      
      
      ord = match(bgen_order$ID_1[-1], my_data_BMI2$ID)
      colname = paste0("BMI2", overlap)
      Phenofile %>%
        mutate({{colname}} := replace_na(Pheno_res_BMI2[ord], -999)) -> Phenofile
    }
    
    
    
    Phenofile %>%
      setnames(paste0(colnames(.), "_", ID)) -> Phenofile
    return(Phenofile)
    
    
  }
  
  
  my_phenofile <- c(simulation_IDs) %>%
    map(getPhenoFile_OneRepetition) %>%
    reduce(cbind)
  
  write.table(my_phenofile, paste0(paste0(simulation_IDs, collapse="-"),"_phenofile"), sep=" ", quote=F, row.names=F)
  
  
  
  if(JURA){
    my_slurm_options = list(partition = "normal", 
                            time="2-00:00:00",  
                            `cpus-per-task`=8, mem="16G")
  } else {
    my_slurm_options = list(partition = "cluster", 
                            time="2-00:00:00",  
                            `cpus-per-task`=8)
  }
  
  name = paste0(paste0(simulation_IDs, collapse="-"), "_")
  
  pars <- data.frame(chr = c(1:22),
                     phenofile = paste0(paste0(simulation_IDs, collapse="-"),"_phenofile"),
                     folder = getwd(),
                     name,
                     JURA)
  
  # launch it
  sjobGWAS <- slurm_apply(launch_bgenie, pars,
                          jobname = paste0('Run_GWAS_BMI-BMI2-', paste0(simulation_IDs, collapse="_")),
                          nodes = 22, cpus_per_node = 1,
                          # specify where packages are stored
                          libPaths=.libPaths(),
                          # specify partition (sgg/cluster/cluster2)
                          slurm_options = my_slurm_options,
                          submit = TRUE)
  
  
  # Once everything lauched
  # check when it's done and merge results
  wait=T
  while(wait){
    # try once every minute, that's enough
    Sys.sleep(60)
    if(get_job_status(sjobGWAS) == "Job completed or stopped.") wait = F
  }
  
  # read all chromosome files -> deal with _IDs in name
  GWAS <- paste0(name, "chr", 1:22, ".out.gz") %>%
    map(data.table::fread) %>%   
    reduce(rbind)  # 1,151,711    
  
  
  cnames = c("chr", "rsid", "pos", "ref", "alt", "af", "info",
             "beta", "se", "z", "minuslog10p")
  
  
  
  # remove low quality, low AF
  GWAS %>%
    filter(info>0.9, # 1,151,150
           a_0 %in% c("A", "C", "G", "T"),
           a_1 %in% c("A", "C", "G", "T"), # 1,151,150
           pmin(af, 1-af)>0.01) -> GWAS_cleaned # 1,151,058
  
  
  # clean files 
  for(ID in simulation_IDs){
    # BMI
    GWAS_cleaned %>%
      select(1:7, starts_with(paste0("BMI_", ID))) %>%
      setNames(cnames) %>%
      mutate(p=10^-minuslog10p) -> GWAS_BMI
    data.table::fwrite(GWAS_BMI, paste0("GWAS_BMI_", ID, ".tsv"), sep="\t")
    # BMI2
    for(overlap in c(100,95,75,50,25,5,0)){
      GWAS_cleaned %>%
        select(1:7, starts_with(paste0("BMI2", overlap, "_", ID))) %>%
        setNames(cnames) %>%
        mutate(p=10^-minuslog10p) -> GWAS_BMI2
      data.table::fwrite(GWAS_BMI2, paste0("GWAS_BMI2", overlap, "_", ID, ".tsv"), sep="\t")
    }
  }
  
  # clean RDS files
  #system("rm results_*")# done once all jobs for a scenario are done
  
  
}

launch_simulation  <- function(){
  library(data.table)
  library(tidyverse)
  options(datatable.freadf.datatable=FALSE)
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
    setwd("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/BMI-BMI2")
  } else {
    setwd("/data/sgg3/ninon/projects/SampleOverlap/Data/BMI-BMI2")
  }
  
 
  if(JURA){
    my_slurm_options = list(partition = "normal", 
                            time="2-00:00:00",  
                            `cpus-per-task`=1, mem="10G")
  } else {
    my_slurm_options = list(partition = sample(c("cluster", "sgg", "cluster2"), 1), 
                            time="1-00:00:00",  mem="10G",
                            `cpus-per-task`=1)
  }
  # check: should not be any results for these IDs already (use ID_seed.txt)
  if(any(file.exists(paste0(1:100, "_seed.txt")))) stop("Some IDs already used for this scenario.")
  #source("../../../Scripts/simulate_data.R")
  
  
  
  # launch simulations 1:100
  id_tolaunch = 1
  launched = 0
  
  while(id_tolaunch<48){
    
    ids = c(id_tolaunch:min(id_tolaunch+5,100))
    id_tolaunch = max(ids) + 1
    
    
    # params for rslurm (add IDs)
    params = list(simulation_IDs = ids)
    
    
    slurm_call(simulate_data, params,
               jobname = paste0("BMI-BMI2_", launched),
               # specify where packages are stored
               libPaths=.libPaths(),
               # specify partition (sgg/cluster/cluster2)
               slurm_options = my_slurm_options,
               submit = TRUE)
    
    
    
    launched = launched+1
    
    # launch bunch of 3*6, then wait 2 hours
    # if(launched %% 3 == 0)   Sys.sleep(2*60*60)
    Sys.sleep(90*60)
    
  }
  
  # JURA
  id_tolaunch = 25
  launched = 4
  
  while(id_tolaunch<48){
    
    ids = c(id_tolaunch:min(id_tolaunch+5,100))
    id_tolaunch = max(ids) + 1
    
    
    # params for rslurm (add IDs)
    params = list(simulation_IDs = ids)
    
    
    slurm_call(simulate_data, params,
               jobname = paste0("BMI-BMI2_", launched),
               # specify where packages are stored
               libPaths=.libPaths(),
               # specify partition (sgg/cluster/cluster2)
               slurm_options = my_slurm_options,
               submit = TRUE)
    
    
    
    launched = launched+1
    
    # launch bunch of 3*6, then wait 2 hours
    # if(launched %% 3 == 0)   Sys.sleep(2*60*60)
    #Sys.sleep(90*60)
    
  }
  
  
  # clean
  # should wait until all jobs are done!!!
  all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0("BMIBMI_", 1:launched, collapse=",")), intern=TRUE))
  while(length(all_myjobs)>1){
    # test every 5 minutes
    Sys.sleep(5*60)
    all_myjobs <- suppressWarnings(system(paste("squeue -n ", paste0("BMIBMI_", 1:launched, collapse=",")), intern=TRUE))
  }

  # remove all BGENIE outputs
  system("rm *.out.gz")

  # remove all _rslurm folders
  system("rm -rf _rslurm*")
  system("rm results_*.RDS")
  
  # that means that in the folders, only these files are left (for each ID):
  # ID_seed.txt                   
  # GWAS_BMI_ID.tsv
  # GWAS_BMI20_ID.tsv
  # ...
  # GWAS_BMI2100_ID.tsv
  # + phenofile (one for multiple IDs)
  
  
}
#launch_simulation()
