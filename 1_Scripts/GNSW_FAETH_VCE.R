###################################################################################
#
# R-Script to run the pipeline for calculating FAETH scores for PIGS
#
# 
#
#
#
#
#
#
#
#
#
##################################################################################
#clean environment
rm(list=ls())

#upload packages
library(dplyr)
library(tidyr)
library(ggplot2)
###################################################################################
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# define the base path
base.path <- args[1]
# Load Base Functions
source(file.path(base.path,"/1_Scripts/GNSW_FAETH_Functions.R"))

###################################################################################
pheno_dir <- paste0(base.path,"/2_Data")
grm_dir <- paste0(base.path,"/3_GRM")

# detecting FUNCTIONAL ANNOTATIONS
# get the names of all files in the directory
FA_list <- list.files(path = grm_dir,
                    include.dirs = TRUE)

cat("\n############################################################ \n")
cat("############ FUNCTIONAL ANNOTATION LAYERS USED ############# \n")
cat("############################################################ \n")
for(map in FA_list){
   cat(paste0("    >  ", map), "\n")
}
cat("\n")

####################################################################################
#
# Run MTG2 to estimate variance components for Target + HD GRM across traits
#
####################################################################################
# GRM files must be built prior to running this pipeline. 
# Please follow the user guide instructions for folder structure and file naming formats used.
# Details on how to use GCTA (Yang et al. 2011) to build GRM from genotypes (PLINK format) are included in materials.
# Yang et al. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 88(1): 76-82.
####################################################################################
#load the phenotypes file
# decor indicates the phenotype file already decorrelated (cholesky decomposition)
phenotypes <- read.table(file = paste0(pheno_dir,"/pheno.txt"), header = TRUE)
trait_list <- names(phenotypes)[3:ncol(phenotypes)]

#define the folder where MTG2 analyses will run
mtg2_dir <- paste0(base.path,"/4_VCE/")

  for(FA in FA_list){
    #define the path to the functional annotation folder within MTG2 folder
    FA_dir <- paste0(mtg2_dir, "/", FA)
    
    dir.create(FA_dir)
    
    setwd(FA_dir)
    
    for(trait in trait_list){
      
      TRT_dir <- paste0(FA_dir, "/", trait)
      
      dir.create(TRT_dir)
      
      setwd(TRT_dir)
      
      if(!file.exists(file.path(TRT_dir, "complete.check"))) {
        # create and write the phenotype file
        pheno_subset <- phenotypes[, c("ANIMAL", "FAM", trait)]
        write.table(pheno_subset, file = "pheno.txt" ,sep = " ", quote = F, row.names = F, col.names = F)
        
        # copy .fam file to be used by MTG2 as guide for GRM numerical ids
        system(paste0("cp ", pheno_dir, "/pheno.fam", " ", TRT_dir,"/pheno.fam"))
        
        # link GRM files to the MTG2 folder (MTG2 will run on links to the real GRM txt file)
        if(FA != "HD"){
          system(paste0("ln -s ", grm_dir, "/", FA, "/G.grm.bin lnk_G_Target.grm"))
          system(paste0("ln -s ", grm_dir, "/HD/G.grm.bin lnk_G_HD.grm"))
        }else{
          system(paste0("ln -s ", grm_dir, "/HD/G.grm.bin lnk_G_HD.grm"))
        }

        #create the .sh file to run MTG2
        crt_MTG2_bash(bash_name = "run_MTG2.sh", memory = "60000", method = "RUN_MTG2_VC", grm_fmt = "GCTA_mbg", ntask = "4", out_file = FA)

        if(FA != "HD"){
          #create file containing names for GRM files (these are the links from real large GRM txt files)
          crt_matlist(GRM1 = "lnk_G_Target.grm", GRM2 = "lnk_G_HD.grm")
        }else{
          crt_matlist(GRM1 = "lnk_G_HD.grm", GRM2 = "")
        }

        system("sbatch run_MTG2.sh")
     }
   }
  }
