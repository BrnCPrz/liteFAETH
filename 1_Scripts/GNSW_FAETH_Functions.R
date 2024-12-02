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
#    FUNCTIONS SCRIPT
#
##################################################################################


crt_MTG2_bash <- function(bash_name, memory = "20000", ntask = "4", method = "RUN_MTG2_VC", grm_fmt = "GCTA_mbg", out_file = "out") {
  
  if(is.null(bash_name)){ cat( " WARNING:Please provide a bash_name !" )}
  
  sink(file = bash_name, type="output")
  #=============== HEADER SECTION ===============#
  cat("#!/bin/bash\n")
  cat(paste0("#SBATCH --job-name=GNSW_",method, "\n"))
  cat("#SBATCH --time=02:00:00 \n")
  cat("#SBATCH --output=output.txt \n")
  cat("#SBATCH --error=error.txt \n")
  cat(paste0("#SBATCH --mem=", memory,"\n"))
  cat(paste0("#SBATCH --ntasks=", ntask,"\n"))
  cat("#SBATCH --account=HG\n")
  cat("#SBATCH --qos=std\n")
  cat("\n")
  cat("\necho \"Model fitting - GREML (required  file ->  greml.matlist) - ESTIMATE VARIANCE COMPONENTS\" \n")
  if(grm_fmt == "GCTA_mbg"){
    cat(paste0("../../../bin/mtg2 -p pheno.fam -d pheno.txt -mbg greml.matlist -mod 1 -thread ", ntask ," -conv 0.0000001 -out ", out_file, ".out > ", out_file, ".log"))
  }
  if(grm_fmt == "CALC_GRM"){
  cat(paste0("../../../bin/mtg2 -p pheno.fam -d pheno.txt -mg greml.matlist -mod 1 -thread ", ntask ," -conv 0.0000001 -out ", out_file, ".out > ", out_file, ".log"))
  }
  cat("\n")
  cat("wait\n")
  cat("echo $(pwd) > complete.check")
  sink()
}

crt_matlist <- function(GRM1 = "lnk_G_Target.grm", GRM2 = "lnk_G_HD.grm") {
  
  sink(file = "greml.matlist", type="output")
  #=============== HEADER SECTION ===============#
  cat(paste0(GRM1,"\n"))
  cat(paste0(GRM2))
  sink()
}