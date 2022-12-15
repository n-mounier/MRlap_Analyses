#### Check effect of sample overlap on MR estimate ####
################ Using Simulated Data #################

## Note, this script only simulate data
# For each simulation :
# get effectives SNPs and effect sizes for X and Y using 1.5M SNPs
# simulate one sample for X and 7 samples (with various overlaps) for Y (using UKBB genotypes)
# run GWAS (for X / for Y) for these 8 samples

# These functions are directly called when using "launch_simulations(...) to simulate 100 datasets

# Relies on:
# - UKB data (bgen files, sample file, ancestry information to restrict analyses to european ancestry individuals)
# - BGENIE
# - slurm workload managed


#### Main Function - no correlated pleiotropy ####
# parameters :
# pi_x / h2_x -> alpha
# pi_y / h2_y -> gamma
# kappa_x
# kappa_u
# alpha
# X = G * alpha + kappa_x * U + eps_x
# Y = G * gamma + alpha * X + kappa_y * U + eps_y

simulate_data <- function(simulation_IDs,
                          scenario,
                          n_A=20000,
                          n_B=20000,
                          pi_x=0.001, # propotion of SNPs affecting X
                          h2_x=0.4,  # heritability of X
                          pi_y=0.01, # proportion of SNPs affecting Y
                          h2_y=0.2 , # heritability of Y
                          kappa_x=0.3, # effect of U on X
                          kappa_y=0.5, # effect of U in Y
                          alpha=0.2){ # causal effect of X on Y

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


  #### check parameters ####
  # check: no more than 6 simulations at one time
  if(length(simulation_IDs)>6) stop("Maximum number of simultaneous simulations (6) exceeded.")

  # check: maximum number of individuals
  if(n_A + n_B > 370000) stop("Maximum number of individuals (370,000) exceeded.")
  # check: simulations paramaters
  if(h2_x+kappa_x^2>=1) stop("X can not have a variance of 1 if h2_x + kappa_x^2 >= 1.")
  if(h2_y+alpha^2+kappa_y^2+2*alpha*kappa_x*kappa_y>=1) stop("Y can not have a variance of 1 if h2_y + alpha^2 + kappa_y^2 + 2 * alpha * kappa_y * kappa_x >= 1.")

  if(JURA){
    setwd("~/scratch_kuta/nmounier/projects/SampleOverlap/Data/Simulations")
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
  if(any(file.exists(paste0(simulation_IDs, "_seed.txt")))) stop("Some IDs already used for this scenario.")


  #### define functions ####
  # for slurm jobs
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

  # to run GWASs
  launch_bgenie <- function(chr, phenofile, folder, name="", JURA=FALSE){
    cat(folder)
    setwd(folder)

    if(JURA){
      system(paste0("~/data_kuta/bin/bgenie_v1.3/bgenie_v1.3_static1 ",
                    "--bgen ~/scratch_kuta/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", chr, ".bgen ",
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

  # data common to all simulations, only read it once
  # all individuals
  # EA_Covariates : to exclude participants with non european ancestry
  # withdrawn : to exclude withdrawn participants
  if(JURA){
    sample_file = data.table::fread("~/data_kuta/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
    EA_Covariates = data.table::fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/SelectedIndividuals_Covariates")
    withdrawn = data.table::fread("~/data_kuta/w19655_20200204.csv")
  } else {
    sample_file = data.table::fread("/data/sgg3/data/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
    EA_Covariates = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Data/Covariates/SelectedIndividuals_Covariates")
    withdrawn = data.table::fread("/data/sgg3/data/UKBB/w19655_20200204.csv")
  }
  sample_file %>%
    slice(-1) -> bgen_order
  bgen_order %>%
    filter(ID_1 %in% EA_Covariates[,1],
           !ID_1 %in% withdrawn[,1]) %>%
    pull(ID_1) -> IDs

  # get SNP info
  rsids = data.table::fread("../../bgen_subsets/subset_SNPs.tsv")
  ## if(JURA){
  ##   rsids = data.table::fread("~/scratch_kuta/nmounier/projects/Run_GWASs/Results/PhysicalActivity/GWAS_Acc425.tsv", select=1:5)
  ## } else{
  ##   rsids = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Results/PhysicalActivity/GWAS_Acc425.tsv", select=1:5)
  ## }

  # to simulate X and Y for one simulation
  simulateData_OneRepetition <- function(ID,
                                         scenario,
                                         n_A=20000,
                                         n_B=20000,
                                         pi_x=0.001, # propotion of SNPs affecting X
                                         h2_x=0.4,  # heritability of X
                                         pi_y=0.01, # propotion of SNPs affecting Y
                                         h2_y=0.2,  # heritability of Y
                                         kappa_x=0.3, # effect of U on X
                                         kappa_y=0.5, # effect of U in Y
                                         alpha=0.2, # causal effect of X on Y
                                         JURA=FALSE,
                                         res_folder){


    # for reproducibility, keep seed info
    seed = round(runif(1, 1, 10^5))
    print(ID)
    print(paste0("seed: ", seed, "\n"))
    set.seed(seed)
    write.table(seed, paste0(ID, "_seed.txt"), sep=",", quote=F, row.names=F, col.names = F)



    #### 1) get samples ####
    # N individuals for trait X
    IDs_X = sample(IDs, n_A)
    # N invididuals for trait Y
    IDs_Y = sample(IDs[!IDs %in% IDs_X], n_B)
    # for Y, get different overlaps
    # 100, 95, 75, 50, 25, 5, 0
    # deal with n_A > n_B vs n_B > n_A (or n_A == n_B)
    if(n_A >= n_B){
      # here, percentage of overlap = percentage of samples of Y that are also in X
      # note, if n_A == n_B, particular case,  percentage of samples of Y than are also in X
      #                           is equalt to percentage of samples of X than are also in Y
      IDs_100 = sample(IDs_X, n_B)
      IDs_95 = c(sample(IDs_X, n_B*0.95), sample(IDs_Y, n_B*0.05))
      IDs_75 = c(sample(IDs_X, n_B*0.75), sample(IDs_Y, n_B*0.25))
      IDs_50 = c(sample(IDs_X, n_B*0.50), sample(IDs_Y, n_B*0.50))
      IDs_25 = c(sample(IDs_X, n_B*0.25), sample(IDs_Y, n_B*0.75))
      IDs_5 =  c(sample(IDs_X, n_B*0.05), sample(IDs_Y, n_B*0.95))
      IDs_0 =  IDs_Y
    } else {
      # here, percentage of overlap = percentage of samples of X that are also in Y
      IDs_100 = c(IDs_X, sample(IDs_Y, (n_B-n_A)))
      IDs_95 = c(sample(IDs_X, n_A*0.95), sample(IDs_Y, n_B-n_A*0.95))
      IDs_75 = c(sample(IDs_X, n_A*0.75), sample(IDs_Y, n_B-n_A*0.75))
      IDs_50 = c(sample(IDs_X, n_A*0.50), sample(IDs_Y, n_B-n_A*0.50))
      IDs_25 = c(sample(IDs_X, n_A*0.25), sample(IDs_Y, n_B-n_A*0.25))
      IDs_5 =  c(sample(IDs_X, n_A*0.05), sample(IDs_Y, n_B-n_A*0.05))
      IDs_0 =  IDs_Y
    }


    all_IDs = c(IDs_X, IDs_Y)
    # we'll need to know the indices of this individual to extract
    # genotypes from bgen files
    all_IDs_inbgen = match(all_IDs, bgen_order$ID_1)

    #### 2) simulate data ####
    # simulate per SNP effect for X & Y
    N_SNPs = nrow(rsids)
    N_effectiveSNPs_X = round(N_SNPs * pi_x)
    effectiveSNPs_X = data.frame(
      rsid = sample(rsids$rsid, size=N_effectiveSNPs_X, replace=F ),
      effect = rnorm(N_effectiveSNPs_X, 0, sqrt(h2_x/N_effectiveSNPs_X)))
    rsids %>%
      mutate(effect_X = case_when(
        rsid %in% effectiveSNPs_X$rsid ~ effectiveSNPs_X$effect[match(rsid, effectiveSNPs_X$rsid)],
        TRUE ~ 0
      )) -> rsids
    write.table(effectiveSNPs_X, paste0(ID, "_effectiveSNPs_X.csv"), sep=",", quote=F, row.names=F)


    N_effectiveSNPs_Y = round(N_SNPs * pi_y)
    effectiveSNPs_Y = data.frame(
      rsid = sample(rsids$rsid, size=N_effectiveSNPs_Y, replace=F ),
      effect = rnorm(N_effectiveSNPs_Y, 0, sqrt(h2_y/N_effectiveSNPs_Y)))
    rsids %>%
      mutate(effect_Y = case_when(
        rsid %in% effectiveSNPs_Y$rsid ~ effectiveSNPs_Y$effect[match(rsid, effectiveSNPs_Y$rsid)],
        TRUE ~ 0
      )) -> rsids
    write.table(effectiveSNPs_Y, paste0(ID, "_effectiveSNPs_Y.csv"), sep=",", quote=F, row.names=F)


    get_SNPeffect <- function(snp, pos, chr, effect, JURA){
      options("stringsAsFactors" = FALSE)

      print(snp)
      data.frame(chromosome=case_when(
        chr %in% 1:9 ~ paste0("0", chr),
        TRUE ~ as.character(chr)),
        start = pos,
        end=pos) -> ranges
      if(JURA){ # rbgen does not deal with ~/ ... ?
        bgen = rbgen::bgen.load(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", as.numeric(ranges$chromosome), ".bgen"), ranges)
      }else{
        bgen = rbgen::bgen.load(paste0("/data/sgg3/ninon/projects/SampleOverlap/Data/bgen_subsets/chr", as.numeric(ranges$chromosome), ".bgen"), ranges)
      }
      # for each individual, probability of g==0, g==1, g==2
      geno = as.data.frame(bgen$data[snp,,])
      geno %>%
        setNames(c("g0", "g1", "g2")) %>%
        rownames_to_column("ID") %>%
        # get genotype from probabilities
        mutate(g = g1*1 + g2*2) %>%
        # standardise genotype
        mutate(gst = scale(g)) -> geno
      # extract relevant individuals
      # here individuals are called "(anonymous_sample_1)" "(anonymous_sample_2)" ...
      # so we need to use "bgen_order" to know which line we'll keep
      x = geno[all_IDs_inbgen,"gst"] * effect
      return(x)
    }

    rsids %>%
      filter(rsid %in% effectiveSNPs_X$rsid) -> pars_X
    pars_X %>%
      transmute(snp =rsid,
                pos=pos,
                chr=chr,
                effect=effect_X,
                JURA=JURA) -> pars_X

    rsids %>%
      filter(rsid %in% effectiveSNPs_Y$rsid) -> pars_Y
    pars_Y %>%
      transmute(snp =rsid,
                pos=pos,
                chr=chr,
                effect=effect_Y,
                JURA=JURA) -> pars_Y

    if(JURA){
      my_slurm_options = list(partition = "normal",
                              time="1-00:00:00",
                              `cpus-per-task`=1,
                              mem="2G")
    } else {
      my_slurm_options = list(partition = "sgg",
                              time="1-00:00:00",
                              `cpus-per-task`=1)
    }

    sjob_getEffects_X <- slurm_apply(get_SNPeffect, pars_X,
                                     jobname = paste0('get_SNPeffect_X_', scenario, "_", ID),
                                     add_objects = c("all_IDs_inbgen"),
                                     nodes = 300, cpus_per_node = 1,
                                     # specify where packages are stored
                                     libPaths=.libPaths(),
                                     # specify partition (sgg/cluster/cluster2)
                                     slurm_options = my_slurm_options,
                                     submit = TRUE)


    if(nrow(pars_Y)>=1){
      sjob_getEffects_Y <- slurm_apply(get_SNPeffect, pars_Y,
                                       jobname = paste0('get_SNPeffect_Y_', scenario, "_", ID),
                                       add_objects = c("all_IDs_inbgen"),
                                       nodes = 100, cpus_per_node = 1,
                                       # specify where packages are stored
                                       libPaths=.libPaths(),
                                       # specify partition (sgg/cluster/cluster2)
                                       slurm_options = my_slurm_options,
                                       submit = TRUE)
    }

    wait=T
    while(wait){
      # try once every minute, that's enough
      Sys.sleep(60)
      #print(get_job_status(sjob_getEffects))
      if(get_job_status(sjob_getEffects_X) == "Job completed or stopped.") wait = F
    }

    if(nrow(pars_Y)>=1){
      wait=T
      while(wait){
        # try once every minute, that's enough
        Sys.sleep(60)
        #print(get_job_status(sjob_getEffects))
        if(get_job_status(sjob_getEffects_Y) == "Job completed or stopped.") wait = F
      }
    }

    # using rslurm to get the results (i.e. loading all results files
    # into a single objet) uses too much memory
    # X = rslurm::get_slurm_out(sjob_getEffects, "table")

    # -> load each object at once, and get the sum (only 1 value / individual)
    # before loading the next one
    setwd(paste0('_rslurm_', sjob_getEffects_X$jobname))
    res_files = list.files(pattern = "results*")
    N=n_A+n_B
    X = rep(0, N)
    for(file in res_files){
      slurm_out = readRDS(file)
      slurm_outd = as.data.frame(slurm_out)
      # slurm_outd -> rows: individuals / columns: SNPs
      X = X + rowSums(slurm_outd)
    }
    setwd("..")

    N=n_A+n_B
    Y = rep(0, N)
    if(nrow(pars_Y)>=1){
      setwd(paste0('_rslurm_', sjob_getEffects_Y$jobname))
      res_files = list.files(pattern = "results*")
      for(file in res_files){
        slurm_out = readRDS(file)
        slurm_outd = as.data.frame(slurm_out)
        # slurm_outd -> rows: individuals / columns: SNPs
        Y = Y + rowSums(slurm_outd)
      }
      setwd("..")
    }
    # check: var(X) should be h2_x at this point
    #        and var(Y) should be h2_y
    # NO -> only if g standardised!!!
    # this is what is done now, so ok
    var(X)
    var(Y)



    # add U and eps_x
    U = rnorm(N, 0, 1)
    X = X + kappa_x * U
    # check: var(X) should be h2_x + kappa_x^2 at this point
    var(X)
    h2_x + kappa_x^2

    # to have var(X)=1, we need var(eps_x) = 1 - h2_x + kappa_u^2
    eps_x = rnorm(N, 0, sqrt(1-(h2_x+kappa_x^2)))
    tau_x = kappa_x * U + eps_x
    X = X + eps_x
    # check: var(X) should be 1
    var(X)


    # calculate Y
    Y = Y +  alpha * X + kappa_y * U

    # check: var(Y) should be h2_y + alpha^2 * var(X) + kappa_y^2 * var(U) + 2*alpha*kappa_y*cov(X,U)
    # cov(X,U) = kappa_x
    # alpha^2 + kappa_y^2 at this point
    var(Y)
    h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y *kappa_x


    eps_y = rnorm(N, 0, sqrt(1- (h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y *kappa_x)))
    tau_y = kappa_y * U + eps_y
    Y = Y + eps_y
    # check: var(Y) should be 1
    var(Y)


    # here, also save cor(X,Y) + expected value
    #write.table(cor(X,Y), paste0(res_folder, "/", ID, "_PhenotypicCorrelation.txt"), sep=",", quote=F, row.names=F, col.names = F)

    # create pheno file
    # X : X is pheno for all_IDs
    # order it for IDs_X
    X_pheno = X[match(IDs_X, all_IDs)]
    ord_X = match(bgen_order$ID_1, IDs_X)

    #table(bgen_order$ID_1 == IDs_X[ord_X])

    bgen_order%>%
      transmute(ID=ID_1,
                X = replace_na(X_pheno[ord_X], -999)) %>%
      mutate(ID=NULL)  -> Phenofile


    # Y
    for(overlap in c(100,95,75,50,25,5,0)){

      my_IDs = get(paste0("IDs_", overlap))
      # Y : Y is pheno for all_IDs
      # order it for my_IDs
      Y_pheno = Y[match(my_IDs, all_IDs)]
      ord = match(bgen_order$ID_1,my_IDs)

      colname = paste0("Y", overlap)
      Phenofile %>%
        mutate({{colname}} := replace_na(Y_pheno[ord], -999)) -> Phenofile
    }



    Phenofile %>%
      setnames(paste0(colnames(.), "_", ID)) -> Phenofile
    return(Phenofile)


  }

  #### simulate data for all IDs ####
  simulation_IDs %>%
    map(~simulateData_OneRepetition(., scenario, n_A, n_B, pi_x, h2_x, pi_y, h2_y, kappa_x, kappa_y, alpha, JURA, res_folder)) %>%
    reduce(cbind) -> my_phenofile

  write.table(my_phenofile, paste0(paste0(simulation_IDs, collapse="-"),"_phenofile"), sep=" ", quote=F, row.names=F)


  if(JURA){
    my_slurm_options = list(partition = "normal",
                            time="2-00:00:00",
                            `cpus-per-task`=8,
                            mem="16G")
  } else {
    my_slurm_options = list(partition = "sgg",
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
                          jobname = paste0('Run_GWAS_XY_', scenario, "-", paste0(simulation_IDs, collapse="_")),
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
    # X
    GWAS_cleaned %>%
      select(1:7, starts_with(paste0("X_", ID))) %>%
      setNames(cnames) %>%
      mutate(p=10^-minuslog10p) -> GWAS_X
    data.table::fwrite(GWAS_X, paste0("GWAS_X_", ID, ".tsv"), sep="\t")
    # Y
    for(overlap in c(100,95,75,50,25,5,0)){
      GWAS_cleaned %>%
        select(1:7, starts_with(paste0("Y", overlap, "_", ID))) %>%
        setNames(cnames) %>%
        mutate(p=10^-minuslog10p) -> GWAS_Y
      data.table::fwrite(GWAS_Y, paste0("GWAS_Y", overlap, "_", ID, ".tsv"), sep="\t")
    }
  }

  # clean RDS files + GWAS Results
  #system("rm results_*")
  #system(paste0("rm ", name, "chr*")) # done once all jobs for a scenario are done

}




#### Main Function - correlated pleiotropy ####
# parameters :
# pi_x / h2_x -> gamma_x
# pi_y / h2_y -> gamma_y
# kappa_x
# kappa_u
# alpha
# pi_u / h2_u -> U_g
# q_x / q_y
# X = G * gamma_x + kappa_x * U_e + q_x * U_g + eps_x
# Y = G * gamma_y + alpha * X + kappa_y * U_e + q_y * U_g + eps_y

simulate_data_genU <- function(simulation_IDs,
                               scenario,
                               n_A=20000,
                               n_B=20000,
                               pi_x=0.001, # propotion of SNPs affecting X
                               h2_x=0.4,  # heritability of X
                               pi_y=0.002, # proportion of SNPs affecting Y
                               h2_y=0.2 , # heritability of Y
                               pi_xy=0, # proportion of X SNPs also affecting Y, uncorrelated pleiotropy
                               kappa_x=0.3, # effect of U_e on X
                               kappa_y=0.5, # effect of U_e in Y
                               pi_u=0.005, # propotion of SNPs affecting U_g
                               h2_u=0.1,  # heritability of U_g
                               q_x = 0.1, # effect of U_g on X
                               q_y = 0.2, # effect of U_g on Y
                               alpha=0.2){ # causal effect of X on Y

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


  #### check parameters ####
  # check: no more than 6 simulations at one time
  if(length(simulation_IDs)>6) stop("Maximum number of simultaneous simulations (6) exceeded.")

  # check: maximum number of individuals
  if(n_A + n_B > 370000) stop("Maximum number of individuals (370,000) exceeded.")
  # check: simulations paramaters
  if(h2_x+kappa_x^2>=1) stop("X can not have a variance of 1 if h2_x + kappa_x^2 >= 1.")
  if(h2_y+alpha^2+kappa_y^2+2*alpha*kappa_x*kappa_y>=1) stop("Y can not have a variance of 1 if h2_y + alpha^2 + kappa_y^2 + 2 * alpha * kappa_y * kappa_x >= 1.")

  if(JURA){
    setwd("~/scratch_kuta/nmounier/projects/SampleOverlap/Data/Simulations")
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
  if(any(file.exists(paste0(simulation_IDs, "_seed.txt")))) stop("Some IDs already used for this scenario.")


  #### define functions ####
  # for slurm jobs
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

  # to run GWASs
  launch_bgenie <- function(chr, phenofile, folder, name="", JURA=FALSE){
    cat(folder)
    setwd(folder)

    if(JURA){
      system(paste0("~/data_kuta/bin/bgenie_v1.3/bgenie_v1.3_static1 ",
                    "--bgen ~/scratch_kuta/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", chr, ".bgen ",
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

  # data common to all simulations, only read it once
  # all individuals
  # EA_Covariates : to exclude participants with non european ancestry
  # withdrawn : to exclude withdrawn participants
  if(JURA){
    sample_file = data.table::fread("~/data_kuta/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
    EA_Covariates = data.table::fread("~/scratch_kuta/nmounier/projects/Run_GWASs/Data/Covariates/SelectedIndividuals_Covariates")
    withdrawn = data.table::fread("~/data_kuta/w19655_20200204.csv")
  } else {
    sample_file = data.table::fread("/data/sgg3/data/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
    EA_Covariates = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Data/Covariates/SelectedIndividuals_Covariates")
    withdrawn = data.table::fread("/data/sgg3/data/UKBB/w19655_20200204.csv")
  }
  sample_file %>%
    slice(-1) -> bgen_order
  bgen_order %>%
    filter(ID_1 %in% EA_Covariates[,1],
           !ID_1 %in% withdrawn[,1]) %>%
    pull(ID_1) -> IDs

  # get SNP info
  rsids = data.table::fread("../../bgen_subsets/subset_SNPs.tsv")
  ## if(JURA){
  ##   rsids = data.table::fread("~/scratch_kuta/nmounier/projects/Run_GWASs/Results/PhysicalActivity/GWAS_Acc425.tsv", select=1:5)
  ## } else{
  ##   rsids = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Results/PhysicalActivity/GWAS_Acc425.tsv", select=1:5)
  ## }

  # to simulate X and Y for one simulation
  simulateData_OneRepetition <- function(ID,
                                         scenario,
                                         n_A=20000,
                                         n_B=20000,
                                         pi_x=0.001, # propotion of SNPs affecting X
                                         h2_x=0.4,  # heritability of X
                                         pi_y=0.002, # propotion of SNPs affecting Y
                                         h2_y=0.2,  # heritability of Y
                                         pi_xy,
                                         kappa_x=0.3, # effect of U_e on X
                                         kappa_y=0.5, # effect of U_e in Y
                                         pi_u=0.005, # propotion of SNPs affecting U_g
                                         h2_u=0.1,  # heritability of U_g
                                         q_x = 0.1, # effect of U_g on X
                                         q_y = 0.2, # effect of U_g on Y
                                         alpha=0.2, # causal effect of X on Y
                                         JURA=FALSE,
                                         res_folder){


    # for reproducibility, keep seed info
    seed = round(runif(1, 1, 10^5))
    print(ID)
    print(paste0("seed: ", seed, "\n"))
    set.seed(seed)
    write.table(seed, paste0(ID, "_seed.txt"), sep=",", quote=F, row.names=F, col.names = F)



    #### 1) get samples ####
    # N individuals for trait X
    IDs_X = sample(IDs, n_A)
    # N invididuals for trait Y
    IDs_Y = sample(IDs[!IDs %in% IDs_X], n_B)
    # for Y, get different overlaps
    # 100, 95, 75, 50, 25, 5, 0
    # deal with n_A > n_B vs n_B > n_A (or n_A == n_B)
    if(n_A >= n_B){
      # here, percentage of overlap = percentage of samples of Y that are also in X
      # note, if n_A == n_B, particular case,  percentage of samples of Y than are also in X
      #                           is equalt to percentage of samples of X than are also in Y
      IDs_100 = sample(IDs_X, n_B)
      IDs_95 = c(sample(IDs_X, n_B*0.95), sample(IDs_Y, n_B*0.05))
      IDs_75 = c(sample(IDs_X, n_B*0.75), sample(IDs_Y, n_B*0.25))
      IDs_50 = c(sample(IDs_X, n_B*0.50), sample(IDs_Y, n_B*0.50))
      IDs_25 = c(sample(IDs_X, n_B*0.25), sample(IDs_Y, n_B*0.75))
      IDs_5 =  c(sample(IDs_X, n_B*0.05), sample(IDs_Y, n_B*0.95))
      IDs_0 =  IDs_Y
    } else {
      # here, percentage of overlap = percentage of samples of X that are also in Y
      IDs_100 = c(IDs_X, sample(IDs_Y, (n_B-n_A)))
      IDs_95 = c(sample(IDs_X, n_A*0.95), sample(IDs_Y, n_B-n_A*0.95))
      IDs_75 = c(sample(IDs_X, n_A*0.75), sample(IDs_Y, n_B-n_A*0.75))
      IDs_50 = c(sample(IDs_X, n_A*0.50), sample(IDs_Y, n_B-n_A*0.50))
      IDs_25 = c(sample(IDs_X, n_A*0.25), sample(IDs_Y, n_B-n_A*0.25))
      IDs_5 =  c(sample(IDs_X, n_A*0.05), sample(IDs_Y, n_B-n_A*0.05))
      IDs_0 =  IDs_Y
    }


    all_IDs = c(IDs_X, IDs_Y)
    # we'll need to know the indices of this individual to extract
    # genotypes from bgen files
    all_IDs_inbgen = match(all_IDs, bgen_order$ID_1)

    #### 2) simulate data ####
    # simulate per SNP effect for X & Y
    if(pi_xy>0){
      N_SNPs = nrow(rsids)

      N_effectiveSNPs_XY = round(N_SNPs * pi_xy)
      snps_XY = sample(rsids$rsid, size=N_effectiveSNPs_XY, replace=F )
      my_snps = rsids$rsid[! rsids$rsid %in% snps_XY]

      N_effectiveSNPs_Xonly = round(N_SNPs * (pi_x-pi_xy))
      N_effectiveSNPs_X = N_effectiveSNPs_Xonly + N_effectiveSNPs_XY
      effectiveSNPs_X = data.frame(
        rsid = c(sample(my_snps, size=N_effectiveSNPs_Xonly, replace=F ), snps_XY),
        effect = rnorm(N_effectiveSNPs_X, 0, sqrt(h2_x/N_effectiveSNPs_X)))
      # update my_snps to remove X only SNPs
      my_snps = rsids$rsid[! rsids$rsid %in% effectiveSNPs_X$rsid]

      rsids %>%
        mutate(effect_X = case_when(
          rsid %in% effectiveSNPs_X$rsid ~ effectiveSNPs_X$effect[match(rsid, effectiveSNPs_X$rsid)],
          TRUE ~ 0
        )) -> rsids
      write.table(effectiveSNPs_X, paste0(ID, "_effectiveSNPs_X.csv"), sep=",", quote=F, row.names=F)


      N_effectiveSNPs_Yonly = round(N_SNPs * (pi_y-pi_xy))
      N_effectiveSNPs_Y = N_effectiveSNPs_Yonly + N_effectiveSNPs_XY
      effectiveSNPs_Y = data.frame(
        rsid = c(sample(rsids$rsid, size=N_effectiveSNPs_Yonly, replace=F ),snps_XY),
        effect = rnorm(N_effectiveSNPs_Y, 0, sqrt(h2_y/N_effectiveSNPs_Y)))
      rsids %>%
        mutate(effect_Y = case_when(
          rsid %in% effectiveSNPs_Y$rsid ~ effectiveSNPs_Y$effect[match(rsid, effectiveSNPs_Y$rsid)],
          TRUE ~ 0
        )) -> rsids
      write.table(effectiveSNPs_Y, paste0(ID, "_effectiveSNPs_Y.csv"), sep=",", quote=F, row.names=F)



    } else { # otherwise, random overlap between X & Y SNPs allowed
      N_SNPs = nrow(rsids)
      N_effectiveSNPs_X = round(N_SNPs * (pi_x))
      effectiveSNPs_X = data.frame(
        rsid = sample(rsids$rsid, size=N_effectiveSNPs_X, replace=F ),
        effect = rnorm(N_effectiveSNPs_X, 0, sqrt(h2_x/N_effectiveSNPs_X)))
      rsids %>%
        mutate(effect_X = case_when(
          rsid %in% effectiveSNPs_X$rsid ~ effectiveSNPs_X$effect[match(rsid, effectiveSNPs_X$rsid)],
          TRUE ~ 0
        )) -> rsids
      write.table(effectiveSNPs_X, paste0(ID, "_effectiveSNPs_X.csv"), sep=",", quote=F, row.names=F)


      N_effectiveSNPs_Y = round(N_SNPs * pi_y)
      effectiveSNPs_Y = data.frame(
        rsid = sample(rsids$rsid, size=N_effectiveSNPs_Y, replace=F ),
        effect = rnorm(N_effectiveSNPs_Y, 0, sqrt(h2_y/N_effectiveSNPs_Y)))
      rsids %>%
        mutate(effect_Y = case_when(
          rsid %in% effectiveSNPs_Y$rsid ~ effectiveSNPs_Y$effect[match(rsid, effectiveSNPs_Y$rsid)],
          TRUE ~ 0
        )) -> rsids
      write.table(effectiveSNPs_Y, paste0(ID, "_effectiveSNPs_Y.csv"), sep=",", quote=F, row.names=F)

    }

    # for U_g
    # U SNPs should be independent
    unusable_SNPs = c(effectiveSNPs_X$rsid, effectiveSNPs_Y$rsid)
    usable_SNPs = rsids$rsid[!rsids$rsid %in% unusable_SNPs]
    N_effectiveSNPs_U = round(N_SNPs * pi_u)
    effectiveSNPs_U = data.frame(
      rsid = sample(usable_SNPs, size=N_effectiveSNPs_U, replace=F ),
      effect = rnorm(N_effectiveSNPs_U, 0, sqrt(h2_u/N_effectiveSNPs_U)))
    write.table(effectiveSNPs_U, paste0(ID, "_effectiveSNPs_U.csv"), sep=",", quote=F, row.names=F)
    rsids %>%
      mutate(effect_U = case_when(
        rsid %in% effectiveSNPs_U$rsid ~ effectiveSNPs_U$effect[match(rsid, effectiveSNPs_U$rsid)],
        TRUE ~ 0
      )) -> rsids


    get_SNPeffect <- function(snp, pos, chr, effect, JURA){
      options("stringsAsFactors" = FALSE)

      print(snp)
      data.frame(chromosome=case_when(
        chr %in% 1:9 ~ paste0("0", chr),
        TRUE ~ as.character(chr)),
        start = pos,
        end=pos) -> ranges
      if(JURA){ # rbgen does not deal with ~/ ... ?
        bgen = rbgen::bgen.load(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", as.numeric(ranges$chromosome), ".bgen"), ranges)
      }else{
        bgen = rbgen::bgen.load(paste0("/data/sgg3/ninon/projects/SampleOverlap/Data/bgen_subsets/chr", as.numeric(ranges$chromosome), ".bgen"), ranges)
      }
      # for each individual, probability of g==0, g==1, g==2
      geno = as.data.frame(bgen$data[snp,,])
      geno %>%
        setNames(c("g0", "g1", "g2")) %>%
        rownames_to_column("ID") %>%
        # get genotype from probabilities
        mutate(g = g1*1 + g2*2) %>%
        # standardise genotype
        mutate(gst = scale(g)) -> geno
      # extract relevant individuals
      # here individuals are called "(anonymous_sample_1)" "(anonymous_sample_2)" ...
      # so we need to use "bgen_order" to know which line we'll keep
      x = geno[all_IDs_inbgen,"gst"] * effect
      return(x)
    }

    rsids %>%
      filter(rsid %in% effectiveSNPs_X$rsid) -> pars_X
    pars_X %>%
      transmute(snp =rsid,
                pos=pos,
                chr=chr,
                effect=effect_X,
                JURA=JURA) -> pars_X

    rsids %>%
      filter(rsid %in% effectiveSNPs_Y$rsid) -> pars_Y
    pars_Y %>%
      transmute(snp =rsid,
                pos=pos,
                chr=chr,
                effect=effect_Y,
                JURA=JURA) -> pars_Y

    rsids %>%
      filter(rsid %in% effectiveSNPs_U$rsid) -> pars_U
    pars_U %>%
      transmute(snp =rsid,
                pos=pos,
                chr=chr,
                effect=effect_U,
                JURA=JURA) -> pars_U

    if(JURA){
      my_slurm_options = list(partition = "normal",
                              time="0-00:30:00",
                              `cpus-per-task`=1,
                              mem="2G")
    } else {
      my_slurm_options = list(partition = "sgg",
                              time="0-00:30:00",
                              `cpus-per-task`=1)
    }

    sjob_getEffects_X <- slurm_apply(get_SNPeffect, pars_X,
                                     jobname = paste0('get_SNPeffect_X_', scenario, "_", ID),
                                     add_objects = c("all_IDs_inbgen"),
                                     nodes = 100, cpus_per_node = 1,
                                     # specify where packages are stored
                                     libPaths=.libPaths(),
                                     # specify partition (sgg/cluster/cluster2)
                                     slurm_options = my_slurm_options,
                                     submit = TRUE)


    if(nrow(pars_Y)>0){
      sjob_getEffects_Y <- slurm_apply(get_SNPeffect, pars_Y,
                                       jobname = paste0('get_SNPeffect_Y_', scenario, "_", ID),
                                       add_objects = c("all_IDs_inbgen"),
                                       nodes = 100, cpus_per_node = 1,
                                       # specify where packages are stored
                                       libPaths=.libPaths(),
                                       # specify partition (sgg/cluster/cluster2)
                                       slurm_options = my_slurm_options,
                                       submit = TRUE)
    }

    sjob_getEffects_U <- slurm_apply(get_SNPeffect, pars_U,
                                     jobname = paste0('get_SNPeffect_U_', scenario, "_", ID),
                                     add_objects = c("all_IDs_inbgen"),
                                     nodes = 100, cpus_per_node = 1,
                                     # specify where packages are stored
                                     libPaths=.libPaths(),
                                     # specify partition (sgg/cluster/cluster2)
                                     slurm_options = my_slurm_options,
                                     submit = TRUE)

    wait=T
    while(wait){
      # try once every minute, that's enough
      Sys.sleep(60)
      #print(get_job_status(sjob_getEffects))
      if(get_job_status(sjob_getEffects_X) == "Job completed or stopped.") wait = F
    }

    if(nrow(pars_Y)>0){
      wait=T
      while(wait){
        # try once every minute, that's enough
        Sys.sleep(60)
        #print(get_job_status(sjob_getEffects))
        if(get_job_status(sjob_getEffects_Y) == "Job completed or stopped.") wait = F
      }
    }

    wait=T
    while(wait){
      # try once every minute, that's enough
      Sys.sleep(60)
      #print(get_job_status(sjob_getEffects))
      if(get_job_status(sjob_getEffects_U) == "Job completed or stopped.") wait = F
    }


    # using rslurm to get the results (i.e. loading all results files
    # into a single objet) uses too much memory
    # X = rslurm::get_slurm_out(sjob_getEffects, "table")

    # -> load each object at once, and get the sum (only 1 value / individual)
    # before loading the next one
    setwd(paste0('_rslurm_', sjob_getEffects_X$jobname))
    res_files = list.files(pattern = "results*")
    N=n_A+n_B
    X = rep(0, N)
    for(file in res_files){
      slurm_out = readRDS(file)
      slurm_outd = as.data.frame(slurm_out)
      # slurm_outd -> rows: individuals / columns: SNPs
      X = X + rowSums(slurm_outd)
    }
    setwd("..")

    Y = rep(0, N)

    if(nrow(pars_Y)>0){
      setwd(paste0('_rslurm_', sjob_getEffects_Y$jobname))
      res_files = list.files(pattern = "results*")
      N=n_A+n_B
      for(file in res_files){
        slurm_out = readRDS(file)
        slurm_outd = as.data.frame(slurm_out)
        # slurm_outd -> rows: individuals / columns: SNPs
        Y = Y + rowSums(slurm_outd)
      }
      setwd("..")
    }
    # check: var(X) should be h2_x at this point
    #        and var(Y) should be h2_y
    # NO -> only if g standardised!!!
    # this is what is done now, so ok
    var(X)
    var(Y)

    # U_g
    setwd(paste0('_rslurm_', sjob_getEffects_U$jobname))
    res_files = list.files(pattern = "results*")
    N=n_A+n_B
    U_g = rep(0, N)
    for(file in res_files){
      slurm_out = readRDS(file)
      slurm_outd = as.data.frame(slurm_out)
      # slurm_outd -> rows: individuals / columns: SNPs
      U_g = U_g + rowSums(slurm_outd)
    }
    setwd("..")

    var(U_g)
    # first, U should be N(0,1), so add some error term
    eps_U = rnorm(N, 0, sqrt(1-(h2_u)))
    U_g = U_g + eps_U
    var(U_g) # should 1 now!


    # add U_e and eps_x
    U_e = rnorm(N, 0, 1)
    X = X + kappa_x * U_e + q_x * U_g
    # check: var(X) should be h2_x + kappa_x^2  + q_x^2 at this point
    var(X)
    h2_x + kappa_x^2 + q_x^2

    # to have var(X)=1, we need var(eps_x) = 1 - h2_x + kappa_u^2
    eps_x = rnorm(N, 0, sqrt(1-(h2_x+kappa_x^2+q_x^2)))
    X = X + eps_x
    # check: var(X) should be 1
    var(X)


    # calculate Y
    Y = Y +  alpha * X + kappa_y * U_e + q_y * U_g
    # check: var(Y) should be h2_y + alpha^2 * var(X) + kappa_y^2 * var(U_e) + 2*alpha*kappa_y*cov(X,U_e) + q_y^2 * var(U_g) + 2*alpha*q_y*cov(X, U_g)
    # cov(X,U_g) = kappa_x
    # cov(X,U_e) = q_x

    var(Y)
    h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y*kappa_x + q_y^2 + 2*alpha*q_x*q_y

    eps_y = rnorm(N, 0, sqrt(1- (h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y*kappa_x + q_y^2 + 2*alpha*q_x*q_y)))
    Y = Y + eps_y
    # check: var(Y) should be 1
    var(Y)


    # here, also save cor(X,Y) + expected value
    #write.table(cor(X,Y), paste0(res_folder, "/", ID, "_PhenotypicCorrelation.txt"), sep=",", quote=F, row.names=F, col.names = F)

    # create pheno file
    # X : X is pheno for all_IDs
    # order it for IDs_X
    X_pheno = X[match(IDs_X, all_IDs)]
    ord_X = match(bgen_order$ID_1, IDs_X)

    #table(bgen_order$ID_1 == IDs_X[ord_X])

    bgen_order%>%
      transmute(ID=ID_1,
                X = replace_na(X_pheno[ord_X], -999)) %>%
      mutate(ID=NULL)  -> Phenofile


    # Y
    for(overlap in c(100,95,75,50,25,5,0)){

      my_IDs = get(paste0("IDs_", overlap))
      # Y : Y is pheno for all_IDs
      # order it for my_IDs
      Y_pheno = Y[match(my_IDs, all_IDs)]
      ord = match(bgen_order$ID_1,my_IDs)

      colname = paste0("Y", overlap)
      Phenofile %>%
        mutate({{colname}} := replace_na(Y_pheno[ord], -999)) -> Phenofile
    }



    Phenofile %>%
      setnames(paste0(colnames(.), "_", ID)) -> Phenofile
    return(Phenofile)


  }

  #### simulate data for all IDs ####
  simulation_IDs %>%
    map(~simulateData_OneRepetition(., scenario, n_A, n_B, pi_x, h2_x, pi_y, h2_y, pi_xy, kappa_x, kappa_y, pi_u, h2_u, q_x, q_y, alpha, JURA, res_folder)) %>%
    reduce(cbind) -> my_phenofile

  write.table(my_phenofile, paste0(paste0(simulation_IDs, collapse="-"),"_phenofile"), sep=" ", quote=F, row.names=F)


  if(JURA){
    my_slurm_options = list(partition = "normal",
                            time="2-00:00:00",
                            `cpus-per-task`=8,
                            mem="16G")
  } else {
    my_slurm_options = list(partition = sample(c("sgg", "cluster2"), 1),
                            time="0-20:00:00",
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
                          jobname = paste0('Run_GWAS_XY_', scenario, "-", paste0(simulation_IDs, collapse="_")),
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
    # X
    GWAS_cleaned %>%
      select(1:7, starts_with(paste0("X_", ID))) %>%
      setNames(cnames) %>%
      mutate(p=10^-minuslog10p) -> GWAS_X
    data.table::fwrite(GWAS_X, paste0("GWAS_X_", ID, ".tsv"), sep="\t")
    # Y
    for(overlap in c(100,95,75,50,25,5,0)){
      GWAS_cleaned %>%
        select(1:7, starts_with(paste0("Y", overlap, "_", ID))) %>%
        setNames(cnames) %>%
        mutate(p=10^-minuslog10p) -> GWAS_Y
      data.table::fwrite(GWAS_Y, paste0("GWAS_Y", overlap, "_", ID, ".tsv"), sep="\t")
    }
  }

  # clean RDS files + GWAS Results
  #system("rm results_*")
  #system(paste0("rm ", name, "chr*")) # done once all jobs for a scenario are done

}



#### Main Function - case control (analysed as continuous) ####
# parameters :
# pi_x / h2_x -> gamma_x
# pi_y / h2_y -> gamma_y
# kappa_x
# kappa_u
# alpha
# pi_u / h2_u -> U_g
# q_x / q_y
# prevalence
# X_liab = G * gamma_x + kappa_x * U_e + q_x * U_g + eps_x
# X_CC = 0 or 1 (depending on X_liab and on threshold derived from prevalence)
# Y = G * gamma_y + alpha * X_CC + kappa_y * U_e + q_y * U_g + eps_y
simulate_data_CC <- function(simulation_IDs,
                             scenario,
                             n_A=20000,
                             n_B=20000,
                             pi_x=0.001, # propotion of SNPs affecting X
                             h2_x=0.4,  # heritability of X
                             pi_y=0.01, # proportion of SNPs affecting Y
                             h2_y=0.2 , # heritability of Y
                             kappa_x=0.3, # effect of U on X
                             kappa_y=0.5, # effect of U in Y
                             alpha=0.2, # causal effect of X on Y
                             prevalence=0.1){ # prevalence

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


  #### check parameters ####
  # check: no more than 6 simulations at one time
  if(length(simulation_IDs)>10) stop("Maximum number of simultaneous simulations (6) exceeded.")

  # check: maximum number of individuals
  if(n_A + n_B > 370000) stop("Maximum number of individuals (370,000) exceeded.")
  # check: simulations paramaters
  if(h2_x+kappa_x^2>=1) stop("X can not have a variance of 1 if h2_x + kappa_x^2 >= 1.")
  if(h2_y+alpha^2+kappa_y^2+2*alpha*kappa_x*kappa_y>=1) stop("Y can not have a variance of 1 if h2_y + alpha^2 + kappa_y^2 + 2 * alpha * kappa_y * kappa_x >= 1.")

  if(JURA){
    setwd("~/scratch_kuta/nmounier/projects/SampleOverlap/Data/Simulations")
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
  if(any(file.exists(paste0(simulation_IDs, "_seed.txt")))) stop("Some IDs already used for this scenario.")


  #### define functions ####
  # for slurm jobs
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

  # to run GWASs
  launch_bgenie <- function(chr, phenofile, folder, name="", JURA=FALSE){
    cat(folder)
    setwd(folder)

    if(JURA){
      system(paste0("~/data_kuta/bin/bgenie_v1.3/bgenie_v1.3_static1 ",
                    "--bgen ~/scratch_kuta/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", chr, ".bgen ",
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

  # data common to all simulations, only read it once
  # all individuals
  # EA_Covariates : to exclude participants with non european ancestry
  # withdrawn : to exclude withdrawn participants
  if(JURA){
    sample_file = data.table::fread("~/data_kuta/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
    EA_Covariates = data.table::fread("~/scratch_kuta/nmounier/projects/Run_GWASs/Data/Covariates/SelectedIndividuals_Covariates")
    withdrawn = data.table::fread("~/data_kuta/w19655_20200204.csv")
  } else {
    sample_file = data.table::fread("/data/sgg3/data/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample")
    EA_Covariates = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Data/Covariates/SelectedIndividuals_Covariates")
    withdrawn = data.table::fread("/data/sgg3/data/UKBB/w19655_20200204.csv")
  }
  sample_file %>%
    slice(-1) -> bgen_order
  bgen_order %>%
    filter(ID_1 %in% EA_Covariates[,1],
           !ID_1 %in% withdrawn[,1]) %>%
    pull(ID_1) -> IDs

  # get SNP info
  rsids = data.table::fread("../../bgen_subsets/subset_SNPs.tsv")
  ## if(JURA){
  ##   rsids = data.table::fread("~/scratch_kuta/nmounier/projects/Run_GWASs/Results/PhysicalActivity/GWAS_Acc425.tsv", select=1:5)
  ## } else{
  ##   rsids = data.table::fread("/data/sgg2/ninon/projects/Run_GWASs/Results/PhysicalActivity/GWAS_Acc425.tsv", select=1:5)
  ## }

  # to simulate X and Y for one simulation
  simulateData_OneRepetition <- function(ID,
                                         scenario,
                                         n_A=20000,
                                         n_B=20000,
                                         pi_x=0.001, # propotion of SNPs affecting X
                                         h2_x=0.4,  # heritability of X
                                         pi_y=0.01, # propotion of SNPs affecting Y
                                         h2_y=0.2,  # heritability of Y
                                         kappa_x=0.3, # effect of U on X
                                         kappa_y=0.5, # effect of U in Y
                                         alpha=0.2, # causal effect of X on Y
                                         prevalence=0.1, #prevalence
                                         JURA=FALSE,
                                         res_folder){


    # for reproducibility, keep seed info
    seed = round(runif(1, 1, 10^5))
    print(ID)
    print(paste0("seed: ", seed, "\n"))
    set.seed(seed)
    write.table(seed, paste0(ID, "_seed.txt"), sep=",", quote=F, row.names=F, col.names = F)



    #### 1) get samples ####
    # N individuals for trait X
    IDs_X = sample(IDs, n_A)
    # N invididuals for trait Y
    IDs_Y = sample(IDs[!IDs %in% IDs_X], n_B)
    # for Y, get different overlaps
    # 100, 95, 75, 50, 25, 5, 0
    # deal with n_A > n_B vs n_B > n_A (or n_A == n_B)
    if(n_A >= n_B){
      # here, percentage of overlap = percentage of samples of Y that are also in X
      # note, if n_A == n_B, particular case,  percentage of samples of Y than are also in X
      #                           is equalt to percentage of samples of X than are also in Y
      IDs_100 = sample(IDs_X, n_B)
      IDs_95 = c(sample(IDs_X, n_B*0.95), sample(IDs_Y, n_B*0.05))
      IDs_75 = c(sample(IDs_X, n_B*0.75), sample(IDs_Y, n_B*0.25))
      IDs_50 = c(sample(IDs_X, n_B*0.50), sample(IDs_Y, n_B*0.50))
      IDs_25 = c(sample(IDs_X, n_B*0.25), sample(IDs_Y, n_B*0.75))
      IDs_5 =  c(sample(IDs_X, n_B*0.05), sample(IDs_Y, n_B*0.95))
      IDs_0 =  IDs_Y
    } else {
      # here, percentage of overlap = percentage of samples of X that are also in Y
      IDs_100 = c(IDs_X, sample(IDs_Y, (n_B-n_A)))
      IDs_95 = c(sample(IDs_X, n_A*0.95), sample(IDs_Y, n_B-n_A*0.95))
      IDs_75 = c(sample(IDs_X, n_A*0.75), sample(IDs_Y, n_B-n_A*0.75))
      IDs_50 = c(sample(IDs_X, n_A*0.50), sample(IDs_Y, n_B-n_A*0.50))
      IDs_25 = c(sample(IDs_X, n_A*0.25), sample(IDs_Y, n_B-n_A*0.25))
      IDs_5 =  c(sample(IDs_X, n_A*0.05), sample(IDs_Y, n_B-n_A*0.05))
      IDs_0 =  IDs_Y
    }


    all_IDs = c(IDs_X, IDs_Y)
    # we'll need to know the indices of this individual to extract
    # genotypes from bgen files
    all_IDs_inbgen = match(all_IDs, bgen_order$ID_1)

    #### 2) simulate data ####
    # simulate per SNP effect for X & Y
    N_SNPs = nrow(rsids)
    N_effectiveSNPs_X = round(N_SNPs * pi_x)
    effectiveSNPs_X = data.frame(
      rsid = sample(rsids$rsid, size=N_effectiveSNPs_X, replace=F ),
      effect = rnorm(N_effectiveSNPs_X, 0, sqrt(h2_x/N_effectiveSNPs_X)))
    rsids %>%
      mutate(effect_X = case_when(
        rsid %in% effectiveSNPs_X$rsid ~ effectiveSNPs_X$effect[match(rsid, effectiveSNPs_X$rsid)],
        TRUE ~ 0
      )) -> rsids
    write.table(effectiveSNPs_X, paste0(ID, "_effectiveSNPs_X.csv"), sep=",", quote=F, row.names=F)


    N_effectiveSNPs_Y = round(N_SNPs * pi_y)
    effectiveSNPs_Y = data.frame(
      rsid = sample(rsids$rsid, size=N_effectiveSNPs_Y, replace=F ),
      effect = rnorm(N_effectiveSNPs_Y, 0, sqrt(h2_y/N_effectiveSNPs_Y)))
    rsids %>%
      mutate(effect_Y = case_when(
        rsid %in% effectiveSNPs_Y$rsid ~ effectiveSNPs_Y$effect[match(rsid, effectiveSNPs_Y$rsid)],
        TRUE ~ 0
      )) -> rsids
    write.table(effectiveSNPs_Y, paste0(ID, "_effectiveSNPs_Y.csv"), sep=",", quote=F, row.names=F)


    get_SNPeffect <- function(snp, pos, chr, effect, JURA){
      options("stringsAsFactors" = FALSE)

      print(snp)
      data.frame(chromosome=case_when(
        chr %in% 1:9 ~ paste0("0", chr),
        TRUE ~ as.character(chr)),
        start = pos,
        end=pos) -> ranges
      if(JURA){ # rbgen does not deal with ~/ ... ?
        bgen = rbgen::bgen.load(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/nmounier/projects/SampleOverlap/Data/bgen_subsets/chr", as.numeric(ranges$chromosome), ".bgen"), ranges)
      }else{
        bgen = rbgen::bgen.load(paste0("/data/sgg3/ninon/projects/SampleOverlap/Data/bgen_subsets/chr", as.numeric(ranges$chromosome), ".bgen"), ranges)
      }
      # for each individual, probability of g==0, g==1, g==2
      geno = as.data.frame(bgen$data[snp,,])
      geno %>%
        setNames(c("g0", "g1", "g2")) %>%
        rownames_to_column("ID") %>%
        # get genotype from probabilities
        mutate(g = g1*1 + g2*2) %>%
        # standardise genotype
        mutate(gst = scale(g)) -> geno
      # extract relevant individuals
      # here individuals are called "(anonymous_sample_1)" "(anonymous_sample_2)" ...
      # so we need to use "bgen_order" to know which line we'll keep
      x = geno[all_IDs_inbgen,"gst"] * effect
      return(x)
    }

    rsids %>%
      filter(rsid %in% effectiveSNPs_X$rsid) -> pars_X
    pars_X %>%
      transmute(snp =rsid,
                pos=pos,
                chr=chr,
                effect=effect_X,
                JURA=JURA) -> pars_X

    if(nrow(effectiveSNPs_Y)>0){
      rsids %>%
        filter(rsid %in% effectiveSNPs_Y$rsid) -> pars_Y
      pars_Y %>%
        transmute(snp =rsid,
                  pos=pos,
                  chr=chr,
                  effect=effect_Y,
                  JURA=JURA) -> pars_Y
    }

    if(JURA){
      my_slurm_options = list(partition = "normal",
                              time="1-00:00:00",
                              `cpus-per-task`=1,
                              mem="2G")
    } else {
      my_slurm_options = list(partition = "sgg",
                              time="1-00:00:00",
                              `cpus-per-task`=1)
    }

    sjob_getEffects_X <- slurm_apply(get_SNPeffect, pars_X,
                                     jobname = paste0('get_SNPeffect_X_', scenario, "_", ID),
                                     add_objects = c("all_IDs_inbgen"),
                                     nodes = 500, cpus_per_node = 1,
                                     # specify where packages are stored
                                     libPaths=.libPaths(),
                                     # specify partition (sgg/cluster/cluster2)
                                     slurm_options = my_slurm_options,
                                     submit = TRUE)


    if(nrow(effectiveSNPs_Y)>0){

      sjob_getEffects_Y <- slurm_apply(get_SNPeffect, pars_Y,
                                       jobname = paste0('get_SNPeffect_Y_', scenario, "_", ID),
                                       add_objects = c("all_IDs_inbgen"),
                                       nodes = 100, cpus_per_node = 1,
                                       # specify where packages are stored
                                       libPaths=.libPaths(),
                                       # specify partition (sgg/cluster/cluster2)
                                       slurm_options = my_slurm_options,
                                       submit = TRUE)
    }


    wait=T
    while(wait){
      # try once every minute, that's enough
      Sys.sleep(60)
      #print(get_job_status(sjob_getEffects))
      if(get_job_status(sjob_getEffects_X) == "Job completed or stopped.") wait = F
    }

    if(nrow(effectiveSNPs_Y)>0){
      wait=T
      while(wait){
        # try once every minute, that's enough
        Sys.sleep(60)
        #print(get_job_status(sjob_getEffects))
        if(get_job_status(sjob_getEffects_Y) == "Job completed or stopped.") wait = F
      }
    }


    # using rslurm to get the results (i.e. loading all results files
    # into a single objet) uses too much memory
    # X = rslurm::get_slurm_out(sjob_getEffects, "table")

    # -> load each object at once, and get the sum (only 1 value / individual)
    # before loading the next one
    setwd(paste0('_rslurm_', sjob_getEffects_X$jobname))
    res_files = list.files(pattern = "results*")
    N=n_A+n_B
    X = rep(0, N)
    for(file in res_files){
      slurm_out = readRDS(file)
      slurm_outd = as.data.frame(slurm_out)
      # slurm_outd -> rows: individuals / columns: SNPs
      X = X + rowSums(slurm_outd)
    }
    setwd("..")

    N=n_A+n_B
    Y = rep(0, N)
    if(nrow(effectiveSNPs_Y)>0){
      setwd(paste0('_rslurm_', sjob_getEffects_Y$jobname))
      res_files = list.files(pattern = "results*")
      for(file in res_files){
        slurm_out = readRDS(file)
        slurm_outd = as.data.frame(slurm_out)
        # slurm_outd -> rows: individuals / columns: SNPs
        Y = Y + rowSums(slurm_outd)
      }
      setwd("..")
    }
    # check: var(X) should be h2_x at this point
    #        and var(Y) should be h2_y
    # NO -> only if g standardised!!!
    # this is what is done now, so ok
    var(X)
    var(Y)



    # add U and eps_x
    U = rnorm(N, 0, 1)
    X = X + kappa_x * U
    # check: var(X) should be h2_x + kappa_x^2 at this point
    var(X)
    h2_x + kappa_x^2

    # to have var(X)=1, we need var(eps_x) = 1 - h2_x + kappa_u^2
    eps_x = rnorm(N, 0, sqrt(1-(h2_x+kappa_x^2)))
    tau_x = kappa_x * U + eps_x
    X = X + eps_x
    # check: var(X) should be 1
    var(X)

    # then turn X into CC
    threshold = -qnorm(prevalence)
    X_CC = as.numeric(X>threshold)
    var(X_CC)
    mean(X_CC)
    # normalise
    X_CC = scale(X_CC)


    # calculate Y
    Y = Y +  alpha * X_CC + kappa_y * U

    # check: var(Y) should be h2_y + alpha^2 * var(X) + kappa_y^2 * var(U) + 2*alpha*kappa_y*cov(X_CC,U)
    # cov(X,U) = kappa_x
    # but cov(X_CC, U) not anymore
    cov_XCC_U = cov(X_CC, U)
    # cov(X_XX, U) = ???

    # alpha^2 + kappa_y^2 at this point
    var(Y) # 0.3254669
    #h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y * kappa_x
    h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y *  cov_XCC_U


    #eps_y = rnorm(N, 0, sqrt(1- (h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y *kappa_x)))
    eps_y = rnorm(N, 0, sqrt(1- (h2_y + alpha^2 + kappa_y^2 + 2*alpha * kappa_y *cov_XCC_U)))

    Y = Y + eps_y
    # check: var(Y) should be 1
    var(Y)


    # here, also save cor(X,Y) + expected value
    #write.table(cor(X,Y), paste0(res_folder, "/", ID, "_PhenotypicCorrelation.txt"), sep=",", quote=F, row.names=F, col.names = F)

    # create pheno file
    # X : X is pheno for all_IDs
    # order it for IDs_X
    X_pheno = X_CC[match(IDs_X, all_IDs)]
    ord_X = match(bgen_order$ID_1, IDs_X)

    #table(bgen_order$ID_1 == IDs_X[ord_X])

    bgen_order%>%
      transmute(ID=ID_1,
                X = replace_na(X_pheno[ord_X], -999)) %>%
      mutate(ID=NULL)  -> Phenofile


    # Y
    for(overlap in c(100,95,75,50,25,5,0)){

      my_IDs = get(paste0("IDs_", overlap))
      # Y : Y is pheno for all_IDs
      # order it for my_IDs
      Y_pheno = Y[match(my_IDs, all_IDs)]
      ord = match(bgen_order$ID_1,my_IDs)

      colname = paste0("Y", overlap)
      Phenofile %>%
        mutate({{colname}} := replace_na(Y_pheno[ord], -999)) -> Phenofile
    }



    Phenofile %>%
      setnames(paste0(colnames(.), "_", ID)) -> Phenofile
    return(Phenofile)


  }

  #### simulate data for all IDs ####
  simulation_IDs %>%
    map(~simulateData_OneRepetition(., scenario, n_A, n_B, pi_x, h2_x, pi_y, h2_y, kappa_x, kappa_y, alpha, prevalence, JURA, res_folder)) %>%
    reduce(cbind) -> my_phenofile

  write.table(my_phenofile, paste0(paste0(simulation_IDs, collapse="-"),"_phenofile"), sep=" ", quote=F, row.names=F)


  if(JURA){
    my_slurm_options = list(partition = "normal",
                            time="2-00:00:00",
                            `cpus-per-task`=8,
                            mem="16G")
  } else {
    my_slurm_options = list(partition = "sgg",
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
                          jobname = paste0('Run_GWAS_XY_', scenario, "-", paste0(simulation_IDs, collapse="_")),
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
    # X
    GWAS_cleaned %>%
      select(1:7, starts_with(paste0("X_", ID, "_")), starts_with(paste0("X_", ID, "-"))) %>%
      setNames(cnames) %>%
      mutate(p=10^-minuslog10p) -> GWAS_X
    data.table::fwrite(GWAS_X, paste0("GWAS_X_", ID, ".tsv"), sep="\t")
    # Y
    for(overlap in c(100,95,75,50,25,5,0)){
      GWAS_cleaned %>%
        select(1:7, starts_with(paste0("Y", overlap, "_", ID, "_")), starts_with(paste0("Y", overlap, "_", ID, "-"))) %>%
        setNames(cnames) %>%
        mutate(p=10^-minuslog10p) -> GWAS_Y
      data.table::fwrite(GWAS_Y, paste0("GWAS_Y", overlap, "_", ID, ".tsv"), sep="\t")
    }
  }

  # clean RDS files + GWAS Results
  #system("rm results_*")
  #system(paste0("rm ", name, "chr*")) # done once all jobs for a scenario are done

}
