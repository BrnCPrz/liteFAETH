# FAETH Score Calculation Pipeline (Lite Version)

## Overview
The FAETH Score Calculation Pipeline is based on genome-wide restricted maximum likelihood (GREML) analysis using multiple Genomic Relationship Matrices (GRMs). 
These matrices are derived from subsets of sequence variants, partitioned using functional and/or evolutionary annotations. For more details, refer to the methodology described by Xiang et al. (2019) ([bioRxiv](https://www.biorxiv.org/content/10.1101/601658v2)).

## Input Requirements
Ensure all required files are placed in their respective directories. Below are the descriptions of the necessary inputs:

### Phenotype and Animal ID Data
- **Directory**: `2_Data/`
1. **pheno.txt**
   A text file containing phenotypes for FAETH score calculation.
   **Format**:
   ```
   ANIMAL_A | ANIMAL_B | PHENOTYPE 1 | PHENOTYPE 2 | ...
   ```
   - Columns:
     - `ANIMAL_A` and `ANIMAL_B`: Duplicate animal IDs (required for MTG2 software compatibility).
     - Phenotype columns: Observations (missing values should be coded as `NA`).
   - Example file included in the pipeline.

2. **pheno.fam**
   A PLINK `.fam` file for identifying animals in the GRM files.
   **Format**:
   ```
   ANIMAL_ID | ANIMAL_ID | 0 | 0 | 0 | -9
   ```

### Target SNP GRMs and High-Density GRM
- **Directory**: `3_GRM/`
1. **Target GRMs (`G.grm`) and Map Files (`genotypes.bim`)**
   - GRMs must represent subsets of SNPs classified as "relevant" based on functional annotation maps.
   - **HD GRM**:
     - Derived from high-density SNP arrays or randomly selected SNPs (~600K SNPs as in Xiang et al., 2019).
   - **GRM Files**: Standard format from [GCTA software](https://cnsgenomics.com/software/gcta/#Overview).File must be named “G.grm.bin”. (Including the GRM file for HD map)
   - **Map File Format**:
     ```
     CHR | CHR_POS | 0 | POS(bp) | ALLELE_A | ALLELE_B
     ```
     - File must be named “genotypes.bim”

2. **Whole Genome Sequence Map (`WGS.bim`)**
   - **Directory**: `2_Data/`
   - Contains SNP positions from the whole-genome sequence dataset.
   - Must match SNP positions in `genotype.bim` files from the GRM folders.
   - **Format**: Same as `genotypes.bim`.

## Pipeline Execution
The pipeline consists of three R scripts, to be executed in the following order:

1. **`GNSW_FAETH_VCE.R`**
   - Estimates variance components for all GRMs (target + HD GRM) across traits.
   - Outputs written to the `4_VCE/` directory.

2. **`GNSW_FAETH_Extract_h2.R`**
   - Extracts per-SNP heritability from variance component estimation results.
   - Outputs written to the `5_FAETH/` directory.

3. **`GNSW_FAETH_Score_Out.R`**
   - Aggregates per-SNP heritability data to calculate FAETH scores.
   - Outputs final scores to the `5_FAETH/` directory.
   - File `pFAETH_Score.txt` contain FAETH scores calculated for all variants.

## Software Dependencies
- **MTG2**: Used for variance component estimation ([MTG2 Documentation](https://sites.google.com/site/honglee0707/mtg2)).
- **GCTA**: For GRM file generation and processing ([GCTA Overview](https://cnsgenomics.com/software/gcta/#Overview)).

## Key Notes
- Ensure input files are formatted correctly and placed in the appropriate directories before running the pipeline.
- The pipeline may fail if:
  - Chromosome or position information is mismatched between files.
  - Missing or incorrectly formatted input files are provided.
