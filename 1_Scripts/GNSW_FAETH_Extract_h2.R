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

#crate a data_frame to store results
results_h2_all <- tibble()
results_n_snp <- tibble()

for(FA in FA_list){
  #define the path to the functional annotation folder within MTG2 folder
  FA_dir <- paste0(mtg2_dir, "/", FA)
  
  for(trait in trait_list){
    
    TRT_dir <- paste0(FA_dir, "/", trait)
    
    setwd(TRT_dir)
    
    # Read in the file as a character vector
    MTG2_result <- readLines(paste0(TRT_dir,"/", FA, ".out"))
    
    # Find the lines containing "h2" using grep()
    h2_line_indices <- grep("h2", MTG2_result)
    
    # Extract the first number following "h2" in each line using regular expressions and sub()
    h2_numbers <- sapply(h2_line_indices, function(i) {
      as.numeric(sub(".*h2\\s+(\\d+\\.\\d+).*", "\\1", MTG2_result[i]))
    })
    
    target_bim <- paste0(base.path,"/3_GRM/", FA)
    # create the BIM file, sorry don't have the plink output so I just create it
    bim_file <- read.table(paste0(target_bim, "/genotypes.bim"), header = FALSE)
    bim_file1 <- bim_file[,c (1,2,2)]
    bim_file1 <- bim_file1[-1,]
    
    # Read in the bim file and count the number of lines
    num_snps <- nrow(bim_file1)
    
    # Calculate the result
    result <- h2_numbers[1]/num_snps
    # Add and save
    bim_file1[,c(4)] <- result
    
    dir.create(paste0(results_dir, "/", FA, "/", trait))
    write.table(bim_file1, paste0(FA,"_", trait, "_SNP_score.txt"), quote = F, row.names = F, col.names = F)
    
    #Create the dataframe
    FA_nSNP <- data.frame(map = FA, n_SNP = num_snps)
    # store number of SNP per functional annotation layer
    results_n_snp <- results_n_snp %>%
                             bind_rows(., FA_nSNP)

    ###########################################################################################
    
    # extract the lines containing "h2"
    h2_lines <- grep("^\\s*h2", MTG2_result, value = TRUE)
    
    # extract the values in the columns starting with "h2"
    h2_values <- gsub("^\\s*h2\\s+", "", h2_lines)  # remove the "h2" label and leading whitespace
    h2_values <- strsplit(h2_values, "\\s+")  # split each line by whitespace
    h2_values <- unlist(h2_values)  # combine the values into a single vector
    
    # convert the values to numeric
    h2_values <- as.numeric(h2_values)

    #Create the dataframe
    h2.trgt <- data.frame(map = FA,
                             trait = trait,
                             h2 = h2_values[1],
                             se = h2_values[2])
    h2.hd <- data.frame(map = "HD",
                             trait = trait,
                             h2 = h2_values[3],
                             se = h2_values[4])
    
    results_h2_all <- results_h2_all %>%
                             bind_rows(., h2.trgt) %>%
                             bind_rows(., h2.hd)
    
  }
}

###############################################################################################
#recover n_snp for HD panel
HD_bim <- paste0(base.path,"/3_GRM/HD")
# create the BIM file, sorry don't have the plink output so I just create it
bim_file <- read.table(paste0(HD_bim, "/genotypes.bim"), header = FALSE)
bim_file1 <- bim_file[,c (1,2,2)]
bim_file1 <- bim_file1[-1,]
    
# Read in the bim file and count the number of lines
num_snps <- nrow(bim_file1)

#Create the dataframe
FA_nSNP <- data.frame(map = "HD", n_SNP = num_snps)
# store number of SNP per functional annotation layer
results_n_snp <- results_n_snp %>%
                        bind_rows(., FA_nSNP)    
################################################################################################

#merge n_snp to h2 results
results_h2_all$n_snp <- results_n_snp$n_SNP[match(results_h2_all$map, results_n_snp$map)]
#calculate per SNP h2 for each trait/functional annotation
results_h2_all$per_SNP_h2 <- results_h2_all$h2 / results_h2_all$n_snp


#prepare n_snp to be outputed
results_n_snp <- results_n_snp[!duplicated(results_n_snp$map),]



# output results for h2
write.table(results_h2_all, paste0(results_dir,"/results_h2_all.txt"), quote = F, row.names = F, col.names = T)
cat(paste("\n#==============================================================================#\n"))
cat(paste("      Heritabilities for Target_Set and HD panel and \n"))     
cat(paste("      SNP based heritabilities across Function annotations and traits \n"))  
cat(paste0("      printed to:  ", results_dir,"/results_h2_all.txt", "\n"))       
cat(paste("\n#==============================================================================#\n\n"))
# output results for n_snp
write.table(results_n_snp, paste0(results_dir,"/results_n_snp.txt"), quote = F, row.names = F, col.names = T)
cat(paste("\n#==============================================================================#\n"))
cat(paste("      Number of SNP per functional annotation map  \n"))     
cat(paste0("      printed to:  ", results_dir,"/results_n_snp.txt", "\n"))       
cat(paste("\n#==============================================================================#\n"))
