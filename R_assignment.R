library(tidyverse)
library(dplyr)
library(janitor)
library(data.table)
library(stringr)
library(ggpubr)
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
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR")) %>%
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

sort_maize <- join_maize %>%
  arrange(as.numeric(Chromosome), as.numeric(Position))
# sorts by chromosome number and position, in ascending order
view(sort_maize)

chr1_maize <- sort_maize %>%
  filter(Chromosome == 1)
view(chr1_maize)
write.table(chr1_maize, file = "chr1_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr2_maize <- sort_maize %>%
  filter(Chromosome == 2)
view(chr2_maize)
write.table(chr2_maize, file = "chr2_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr3_maize <- sort_maize %>%
  filter(Chromosome == 3)
view(chr3_maize)
write.table(chr3_maize, file = "chr3_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr4_maize <- sort_maize %>%
  filter(Chromosome == 4)
view(chr4_maize)
write.table(chr4_maize, file = "chr4_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr5_maize <- sort_maize %>%
  filter(Chromosome == 5)
view(chr5_maize)
write.table(chr5_maize, file = "chr5_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr6_maize <- sort_maize %>%
  filter(Chromosome == 6)
view(chr6_maize)
write.table(chr6_maize, file = "chr6_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr7_maize <- sort_maize %>%
  filter(Chromosome == 7)
view(chr7_maize)
write.table(chr7_maize, file = "chr7_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr8_maize <- sort_maize %>%
  filter(Chromosome == 8)
view(chr8_maize)
write.table(chr8_maize, file = "chr8_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr9_maize <- sort_maize %>%
  filter(Chromosome == 9)
view(chr9_maize)
write.table(chr9_maize, file = "chr9_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr10_maize <- sort_maize %>%
  filter(Chromosome == 10)
view(chr10_maize)
write.table(chr10_maize, file = "chr10_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

sort_desc_maize <- join_maize %>%
  arrange(as.numeric(Chromosome), desc(as.numeric(Position))) 
# sorts maize file in ascending order by chromosome and descending order by position
view(sort_desc_maize)

maize_dash <- data.frame(lapply(sort_desc_maize, gsub, pattern = "[?]", replacement = "-"))
# replaces ? with -
view(maize_dash)

chr1_maize_dash <- maize_dash %>%
  filter(Chromosome == 1)
view(chr1_maize_dash)
write.table(chr1_maize_dash, file = "chr1_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr2_maize_dash <- maize_dash %>%
  filter(Chromosome == 2)
view(chr2_maize_dash)
write.table(chr2_maize_dash, file = "chr2_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr3_maize_dash <- maize_dash %>%
  filter(Chromosome == 3)
view(chr3_maize_dash)
write.table(chr3_maize_dash, file = "chr3_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr4_maize_dash <- maize_dash %>%
  filter(Chromosome == 4)
view(chr4_maize_dash)
write.table(chr4_maize_dash, file = "chr4_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr5_maize_dash <- maize_dash %>%
  filter(Chromosome == 5)
view(chr5_maize_dash)
write.table(chr5_maize_dash, file = "chr5_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr6_maize_dash <- maize_dash %>%
  filter(Chromosome == 6)
view(chr6_maize_dash)
write.table(chr6_maize_dash, file = "chr6_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr7_maize_dash <- maize_dash %>%
  filter(Chromosome == 7)
view(chr7_maize_dash)
write.table(chr7_maize_dash, file = "chr7_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr8_maize_dash <- maize_dash %>%
  filter(Chromosome == 8)
view(chr8_maize_dash)
write.table(chr8_maize_dash, file = "chr8_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr9_maize_dash <- maize_dash %>%
  filter(Chromosome == 9)
view(chr9_maize_dash)
write.table(chr9_maize_dash, file = "chr9_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr10_maize_dash <- maize_dash %>%
  filter(Chromosome == 10)
view(chr10_maize_dash)
write.table(chr10_maize_dash, file = "chr10_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

##### Data processing: teosinte files
geno_teosinte <- geno %>%
  filter(Group %in% c("ZMPBA", "ZMPIL", "ZMPJA")) %>% 
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
view(geno_teosinte)

join_teosinte <- snp_cut %>%
  inner_join(geno_teosinte, by = "SNP_ID")
# joins both files - inner_join makes sure we take only SNP_IDs in both files
view(join_teosinte)

sort_teosinte <- join_teosinte %>%
  arrange(as.numeric(Chromosome), as.numeric(Position))
# sorts by chromosome number and position, in ascending order
view(sort_teosinte)

chr1_teosinte <- sort_teosinte %>%
  filter(Chromosome == 1)
view(chr1_teosinte)
write.table(chr1_teosinte, file = "chr1_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr2_teosinte <- sort_teosinte %>%
  filter(Chromosome == 2)
view(chr2_teosinte)
write.table(chr2_teosinte, file = "chr2_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr3_teosinte <- sort_teosinte %>%
  filter(Chromosome == 3)
view(chr3_teosinte)
write.table(chr3_teosinte, file = "chr3_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr4_teosinte <- sort_teosinte %>%
  filter(Chromosome == 4)
view(chr4_teosinte)
write.table(chr4_teosinte, file = "chr4_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr5_teosinte <- sort_teosinte %>%
  filter(Chromosome == 5)
view(chr5_teosinte)
write.table(chr5_teosinte, file = "chr5_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr6_teosinte <- sort_teosinte %>%
  filter(Chromosome == 6)
view(chr6_teosinte)
write.table(chr6_teosinte, file = "chr6_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr7_teosinte <- sort_teosinte %>%
  filter(Chromosome == 7)
view(chr7_teosinte)
write.table(chr7_teosinte, file = "chr7_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr8_teosinte <- sort_teosinte %>%
  filter(Chromosome == 8)
view(chr8_teosinte)
write.table(chr8_teosinte, file = "chr8_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr9_teosinte <- sort_teosinte %>%
  filter(Chromosome == 9)
view(chr9_teosinte)
write.table(chr9_teosinte, file = "chr9_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr10_teosinte <- sort_teosinte %>%
  filter(Chromosome == 10)
view(chr10_teosinte)
write.table(chr10_teosinte, file = "chr10_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

sort_desc_teosinte <- join_teosinte %>%
  arrange(as.numeric(Chromosome), desc(as.numeric(Position))) 
# sorts teosinte file in ascending order by chromosome and descending order by position
view(sort_desc_teosinte)

teosinte_dash <- data.frame(lapply(sort_desc_teosinte, gsub, pattern = "[?]", replacement = "-"))
# replaces ? with -
view(teosinte_dash)

chr1_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 1)
view(chr1_teosinte_dash)
write.table(chr1_teosinte_dash, file = "chr1_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr2_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 2)
view(chr2_teosinte_dash)
write.table(chr2_teosinte_dash, file = "chr2_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr3_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 3)
view(chr3_teosinte_dash)
write.table(chr3_teosinte_dash, file = "chr3_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr4_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 4)
view(chr4_teosinte_dash)
write.table(chr4_teosinte_dash, file = "chr4_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr5_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 5)
view(chr5_teosinte_dash)
write.table(chr5_teosinte_dash, file = "chr5_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr6_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 6)
view(chr6_teosinte_dash)
write.table(chr6_teosinte_dash, file = "chr6_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr7_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 7)
view(chr7_teosinte_dash)
write.table(chr7_teosinte_dash, file = "chr7_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr8_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 8)
view(chr8_teosinte_dash)
write.table(chr8_teosinte_dash, file = "chr8_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr9_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 9)
view(chr9_teosinte_dash)
write.table(chr9_teosinte_dash, file = "chr9_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

chr10_teosinte_dash <- teosinte_dash %>%
  filter(Chromosome == 10)
view(chr10_teosinte_dash)
write.table(chr10_teosinte_dash, file = "chr10_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

##### Visualization


sort_maize_bin <- sort_maize %>% 
  filter(Position != "unknown") %>%
  filter(Position != "multiple") %>% # keep only chromosomes of known position
  mutate(position_binned = cut(as.numeric(Position), 10)) %>% # make 10 bins per chromosome
  filter(Chromosome != "multiple") %>% 
  filter(Chromosome != "unknown") %>% # keep only SNPs where chromosome is known
  mutate(Chromosome = fct_relevel(Chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")))
# final step sets levels so graph reads in order of chromosome number

maize <- ggplot(data = sort_maize_bin) + 
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome)) + # colored by chromosome
  xlab("Position(basepairs)") +
  ylab("SNP Distribution") +
  theme(axis.text.x = element_text(angle = 90)) + # changes angle so position is readable
  facet_grid(~as.numeric(Chromosome)) # facet by chromosome


sort_teosinte_bin <- sort_teosinte %>% 
  filter(Position != "unknown") %>%
  filter(Position != "multiple") %>% # keep only chromosomes of known position
  mutate(position_binned = cut(as.numeric(Position), 10)) %>% # make 10 bins per chromosome
  filter(Chromosome != "multiple") %>%
  filter(Chromosome != "unknown") %>% # keep only SNPs where chromosome is known
  mutate(Chromosome = fct_relevel(Chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")))
# final step sets levels so graph reads in order of chromosome number

teosinte <- ggplot(data = sort_teosinte_bin) + 
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome)) + # colored by chromosome
  xlab("Position(basepairs)") +
  ylab("SNP Distribution") +
  theme(axis.text.x = element_text(angle = 90)) + # changes angle so position is readable
  facet_grid(~as.numeric(Chromosome)) # facet by chromosome

SNP_figure <- ggarrange(maize, teosinte, # from ggpubr, arranges graphs
                    labels = c("maize", "teosinte"),
                    ncol = 1, nrow = 2) # graphs areanged with maize on top, teosinte below
SNP_figure

ggsave("SNP_figure.pdf", plot = SNP_figure, width = 12, height = 15, units = "in", dpi = 300) # print to pdf

ggsave("SNP_figure.png", plot = SNP_figure, width = 12, height = 15, units = "in", dpi = 300) # print to png


#SECOND SNP graph starts below
unknown_maize <- sort_maize %>%
  filter(Chromosome == "unknown")
view(unknown_maize)

multiple_maize <- sort_maize %>%
  filter(Chromosome == "multiple")
view(multiple_maize)

mchr1 <- nrow(chr1_maize)
mchr2 <- nrow(chr2_maize)
mchr3 <- nrow(chr3_maize)
mchr4 <- nrow(chr4_maize)
mchr5 <- nrow(chr5_maize)
mchr6 <- nrow(chr6_maize)
mchr7 <- nrow(chr7_maize)
mchr8 <- nrow(chr8_maize)
mchr9 <- nrow(chr9_maize)
mchr10 <- nrow(chr10_maize)
mchrmult <- nrow(multiple_maize)
mchrunkn <- nrow(unknown_maize)

unknown_teosinte <- sort_teosinte %>%
  filter(Chromosome == "unknown")
view(unknown_teosinte)

multiple_teosinte <- sort_teosinte %>%
  filter(Chromosome == "multiple")
view(multiple_teosinte)

tchr1 <- nrow(chr1_teosinte)
tchr2 <- nrow(chr2_teosinte)
tchr3 <- nrow(chr3_teosinte)
tchr4 <- nrow(chr4_teosinte)
tchr5 <- nrow(chr5_teosinte)
tchr6 <- nrow(chr6_teosinte)
tchr7 <- nrow(chr7_teosinte)
tchr8 <- nrow(chr8_teosinte)
tchr9 <- nrow(chr9_teosinte)
tchr10 <- nrow(chr10_teosinte)
tchrmult <- nrow(multiple_teosinte)
tchrunkn <- nrow(unknown_teosinte)

SNP_pos_df <- data.frame(chromosome = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown"),
                    species = c("maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte"),
                    SNP_positions = c(mchr1, mchr2, mchr3, mchr4, mchr5, mchr6, mchr7, mchr8, mchr9, mchr10, mchrmult, mchrunkn, tchr1, tchr2, tchr3, tchr4, tchr5, tchr6, tchr7, tchr8, tchr9, tchr10, tchrmult, tchrunkn))

view(SNP_pos_df)

SNP_pos_df <- SNP_pos_df %>%
  mutate(chromosome = fct_relevel(chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")))

SNP_pos_bar <- ggplot(data = SNP_pos_df, aes(x = species, y = SNP_positions, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Species") +
  ylab("Number of SNP positions") +
  facet_grid(~chromosome)
SNP_pos_bar

ggsave("SNP_pos_bar.pdf", plot = SNP_pos_bar, width = 10, height = 5, units = "in", dpi = 300) # print to pdf

ggsave("SNP_pos_bar.png", plot = SNP_pos_bar, width = 10, height = 5, units = "in", dpi = 300) # print to png


## NEXT graph


geno4 <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t")

geno4$homozygous <- apply(geno4, 1, function(x) length(which(x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G"))))
geno4$heterozygous <- apply(geno4, 1, function(x) length(which(x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("C/G") | x == ("G/A") | x == ("G/C") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G"))))
geno4$missing <- apply(geno4, 1, function(x) length(which(x == ("?/?"))))
geno4$total <- apply(geno4, 1, function(x) length(which(x == ("?/?") | x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G") | x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("C/G") | x == ("G/A") | x == ("G/C") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G"))))

view(geno4)

geno_totals <- geno4 %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "ZMPBA", "ZMPIL", "ZMPJA")) %>%
  mutate(species = fct_recode(Group, 
                             "maize" = "ZMMIL", 
                             "maize" = "ZMMLR",
                             "maize" = "ZMMMR",
                             "teosinte" = "ZMPBA",
                             "teosinte" = "ZMPIL",
                             "teosinte" = "ZMPJA"))
  view(geno_totals)

big_zygo <- geno_totals %>%
  select(Group, homozygous, heterozygous, missing) %>%
  mutate(species = fct_recode(Group, 
                              "maize" = "ZMMIL", 
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA"))


big_zygo_pivot <- big_zygo %>% 
  select(Group, species, homozygous, heterozygous, missing) %>%
  pivot_longer(., cols = c(homozygous, heterozygous, missing), names_to = "zygosity", values_to = "count")
view(big_zygo_pivot)

big_zygo_species_plot <- ggplot(big_zygo_pivot, aes(fill = zygosity, y = count, x = species)) + 
  xlab("Species") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
big_zygo_species_plot


ggsave("zygo_species_plot.pdf", plot = big_zygo_species_plot, width = 5, height = 5, units = "in", dpi = 300) # print to pdf

ggsave("zygo_species_plot.png", plot = big_zygo_species_plot, width = 5, height = 5, units = "in", dpi = 300) # print to png


big_zygo_group_plot <- ggplot(big_zygo_pivot, aes(fill = zygosity, y = count, x = Group)) + 
  xlab("Group") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
big_zygo_group_plot

ggsave("zygo_group_plot.pdf", plot = big_zygo_group_plot, width = 5, height = 5, units = "in", dpi = 300) # print to pdf

ggsave("zygo_group_plot.png", plot = big_zygo_group_plot, width = 5, height = 5, units = "in", dpi = 300) # print to png

## FUN graph

geno5 <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t")

geno5$GC <- apply(geno5, 1, function(x) length(which(x == ("G/C") | x == ("C/G") | x == ("G/G") | x == ("C/C"))))
geno5$not_GC <- apply(geno5, 1, function(x) length(which(x == ("?/?") | x == ("A/A") | x == ("T/T") | x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("G/A") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G"))))
view(geno5)

geno_gc <- geno5 %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "ZMPBA", "ZMPIL", "ZMPJA")) %>%
  mutate(species = fct_recode(Group, 
                              "maize" = "ZMMIL", 
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA"))

gc_pivot <- geno_gc %>% 
  select(Group, species, GC, not_GC) %>%
  pivot_longer(., cols = c(GC, not_GC), names_to = "GC_content", values_to = "count")


view(gc_pivot)

gc_group_plot <- ggplot(gc_pivot, aes(fill = GC_content, y = count, x = Group)) + 
  xlab("Group") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
gc_group_plot

ggsave("gc_group_plot.pdf", plot = gc_group_plot, width = 5, height = 5, units = "in", dpi = 300) # print to pdf

ggsave("gc_group_plot.png", plot = gc_group_plot, width = 5, height = 5, units = "in", dpi = 300) # print to png


gc_species_plot <- ggplot(gc_pivot, aes(fill = GC_content, y = count, x = species)) + 
  xlab("Species") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
gc_species_plot

ggsave("gc_species_plot.pdf", plot = gc_species_plot, width = 5, height = 5, units = "in", dpi = 300) # print to pdf

ggsave("gc_species_plot.png", plot = gc_species_plot, width = 5, height = 5, units = "in", dpi = 300) # print to png
