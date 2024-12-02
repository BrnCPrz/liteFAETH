library(dplyr)

####################################################################################
#
# create genomic relationship matrix for all target_SNP genotypes
#
####################################################################################
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# define the base path
base.path <- args[1]
#define the folder where phenotypes are stored
pheno_dir <- paste0(base.path,"/2_Data")
#define the folder where GRM are stored
grm_dir <- paste0(base.path,"/3_GRM")
#define the folder where MTG2 analyses will run
mtg2_dir <- paste0(base.path,"/4_VCE")
# define folder to store results
results_dir <- paste0(base.path,"/5_FAETH")

# detecting FUNCTIONAL ANNOTATIONS
# get the names of all files in the directory
FA_list <- list.files(path = grm_dir,
                      include.dirs = TRUE)

# decor indicates the phenotype file already decorrelated (cholesky decomposition)
phenotypes <- read.table(file = paste0(pheno_dir,"/pheno.txt"), header = TRUE)
trait_list <- names(phenotypes)[3:ncol(phenotypes)]

cat("\nLOADING WGS map: \n")
WGS.map <- read.table(file = paste0(pheno_dir, "/WGS.bim"), header = F)
names(WGS.map) <- c("CHR", "CHR_POS", "0", "POS", "ALLELE_A", "ALLELE_B")

WGS.map <- WGS.map %>%
            select(CHR, CHR_POS)
cat("FINISHED! \n")
# empty data_frame as container of results
avrg_SNP_h2 <- tibble()

for(FA in FA_list){
  if(FA != "HD"){
  cat(paste("ANNOTATION ->   ", FA, "\n"))
  #define the path to the functional annotation folder within MTG2 folder
  FA_dir <- paste0(mtg2_dir, "/", FA)
  
    for(trait in trait_list){
      
      TRT_dir <- paste0(FA_dir, "/", trait)
      
      setwd(TRT_dir)
      
      # Read in the file as a character vector
      target.map <- read.table(file = paste0(FA, "_", trait, "_SNP_score.txt"), header = F)
      names(target.map) <- c("CHR", "CHR_POS", "CHR_POS_B", "SCORE")
      
      #trim SCORE to only 5 decimals but keep the scientific notation (this saves storage space)
      target.map$SCORE <- as.numeric(format(target.map$SCORE, scientific = TRUE, digits = 5))
      
      # check if SNP-based h2 are negative (artifact from VC estimation)
      neg.flag <- ifelse(mean(target.map$SCORE) < 0, TRUE, FALSE )
      # if SNP-based h2 are negative, correct it to zero
      target.map$SCORE <- ifelse( neg.flag == TRUE, 0, target.map$SCORE)

      target.map <- target.map %>%
                        select(CHR_POS, SCORE)
      
      # Join the two datasets on SNP_ID
      avrg_SNP_h2 <- left_join(WGS.map, target.map, by = "CHR_POS")
      
      #create column name for the scores
      score_ID <- paste0(FA,"_",trait)
      
      # assign scores for that trait to the target SNP in given FA
      WGS.map <- WGS.map %>%
                  mutate(!!score_ID := avrg_SNP_h2$SCORE)
    
  }
 }
}

WGS.map$FAETH_score <- rowMeans(WGS.map[, c(3:ncol(WGS.map))], na.rm = TRUE)

cat("#====================================================================#")
cat(paste("\n\n  Writing pFAETH_Score.bim . . . \n"))
write.table(WGS.map, paste0(results_dir,"/pFAETH_Score.txt"), quote = F, row.names = F, col.names = T)
cat(paste("\n\n  pFAETH_Score.txt  Successfully writen ! \n"))
cat("#====================================================================#")