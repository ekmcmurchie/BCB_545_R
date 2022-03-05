library(tidyverse)
library(dplyr)
library(janitor)
library(data.table)
geno <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t")
snp_pos <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/snp_position.txt", header = TRUE, sep = "\t")
# importing fang_et_al.txt as geno and snp_position.txt as snp_pos

##### Data inspection: fang_et_al.txt
view(geno) # allows us to visually inspect file
# We see that there is a header line with Sample_ID, JG_OTU, and Group
# These three columns are followed by a number of columns labeled with SNP ID
# This header line is followed by the data

ncol(geno) # count number of columns
# 986

dim(geno) # count dimensions (rows and columns)
# 2782 rows (this excludes header - with header, 2783 rows), 986 columns
# Multiplying 2782 by 986 gives the number of "words" excluding header: 2,743,052
# Multiplying 2783 by 986 gives the number of "words" including header: 2,744,038

sum(!is.na(geno)) # this will give word count excluding header
# 2743052

length(geno) # number of columns
# 986

object.size(geno) # size in bytes
# 22681376 bytes

str(geno) # structure of geno
# dataframe, 2782 obs (rows excluding header). of 986 variables (columns)
# Although list output is truncated, first few variables read as characters

##### Data inspection: snp_position.txt
view(snp_pos) # allows us to visually inspect file
# We see that there is a header line with SNP_ID,
# Chromosome, Position, alt_pos, mult_positions, amplicon,
# cdv_map_feature.name, gene,and candidate.random, cvd_marker_id,
# Genaissance_daa_id, Sequenom_daa_id, count_amplicons, count_cmf,and count_gene
# The header is followed by the data

ncol(snp_pos) # count number of columns
# 15 

dim(snp_pos) # count dimensions (rows and columns)
# 983 rows (this excludes header), 15 columns
# Multiplying 983 by 15 gives the number of "words" excluding header: 14,745
# Multiplying 984 by 15 gives the number of "words" including header: 14,760

sum(!is.na(snp_pos)) # this will give word count excluding header
# 14745

length(snp_pos) # number of columns
# 15

object.size(snp_pos) # size in bytes
# 327392 bytes

str(snp_pos) # structure of snp_pos
# dataframe, 983 obs (rows excluding header). of 15 variables (columns)
# SNP_ID, Chromosome, Position, alt_pos, mult_positions, amplicon,
# cdv_map_feature.name, gene,and candidate.random are characters
# cvd_marker_id, Genaissance_daa_id, Sequenom_daa_id, count_amplicons, count_cmf,
# and count_gene are integers

##### Data processing: initial steps
snp_cut <- snp_pos %>%
  select(SNP_ID, Chromosome, Position) # select only SNP_ID, Chromosome, and Position


view(snp_cut)

##### Data processing: maize files
geno_maize <- geno %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMR")) %>%
  # filter by maize groups only, including sample ID
  t() %>%
  # transpose - note it's added automatic column names (V1, V2, etc. in new header)
  as.data.frame() %>%
  # make sure it's a dataframe
  setDT(keep.rownames = TRUE) %>%
  #gives column 1 a placeholder name as well
  row_to_names(row_number = 1) %>%
  # sets first row as column names
  rename(SNP_ID = Sample_ID)
# rename first column as SNP_ID so we can join on like column names
view(geno_maize)

join_maize <- snp_cut %>%
  inner_join(geno_maize, by = "SNP_ID")
# joins both files - inner_join makes sure we take only SNP_IDs in both files
view(join_maize)

