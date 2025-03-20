library(tidyverse)
library(dplyr)
library(readr)
library(gridExtra)

# Load the two files (txt is tab delimited)
snp <- read.csv("00_Data/snp_position.txt", sep = "\t", header = T)
fang <- read.csv("00_Data/fang_et_al_genotypes.txt", sep = "\t", header = T)
fangt <- read.csv("00_Data/transposed_genotypes.txt", sep = "\t", header = T)

# Data Inspection (file size, number of columns, number of lines, etc)
View(fangt) #Look at the files to see if imported correctly
View(snp)
str(fangt) #See the structure of the files
str(snp)
is.data.frame(fangt) #Assignment says import as data frame, double check if they are data frames
is.data.frame(snp)

# Number of rows and columns and size of the files 
nrow(snp)
nrow(fangt)

ncol(snp)
ncol(fangt)

colnames(snp)
colnames(fangt)

format(object.size(snp), units = "auto")
format(object.size(fangt), units = "auto")

# Summary of data inspection 
data_inspection <- data.frame(
  Dataset = c("snp", "fang"),
  Rows = c(nrow(snp), nrow(fang)),
  Columns = c(ncol(snp), ncol(fang)),
  Size_MB = c(
    format(object.size(snp) / 1024^2, digits = 2),
    format(object.size(fangt) / 1024^2, digits = 2)
  )
)

# Data Processing 
# create maize and teosinte subset
maize <- filter(fang, Group %in% c("ZMMIL", "ZMMLR", "ZMMMR"))
maizet <- maize %>% column_to_rownames(., var = "Sample_ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column(., var = "SNP_ID")
maize1 <- maizet[-c(1, 2), ]

teosinte <- filter(fang, Group %in% c("ZMPBA", "ZMPIL", "ZMPJA "))
teosintet <- teosinte %>% column_to_rownames(., var = "Sample_ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column(., var = "SNP_ID")
teosinte1 <- teosintet[-c(1, 2), ]


# select SNP_ID, Chromosome and Position from snp file and merge with maize and teosinte file respectively
join1 <- select(snp, SNP_ID, Chromosome, Position)
joinedmaize <- merge(join1, maize1, by = "SNP_ID", all = TRUE)
joinedteosinte <- merge(join1, teosinte1, by = "SNP_ID", all = TRUE)

# Make directories for the outout files in the output folder 
dir.create("02_Output/./maize_files")
dir.create("02_Output/./teosinte_files")


#encode missing values, replace ? with - but keep both versions (with and without question marks)
joinedmaizerep <- joinedmaize %>% 
 mutate_all(~ gsub("\\?", "-", .))

joinedteosinterep <- joinedteosinte %>% 
  mutate_all(~ gsub("\\?", "-", .))

separate_chr_files <- function(i) {
  maizeinc <- joinedmaize %>% filter(Chromosome == i) %>% arrange(as.numeric(Position))
  maizedec <- joinedmaizerep %>% filter(Chromosome == i) %>% arrange(desc(as.numeric(Position)))
  
  write_tsv(maizeinc, file.path("02_Output/./maize_files/", paste("Maize_chr", i, "increasing.txt", sep = "_")))
  write_tsv(maizedec, file.path("02_Output/./maize_files/", paste("Maize_chr", i, "decreasing.txt", sep = "_")))
  
  teoinc <- joinedteosinte %>% filter(Chromosome == i) %>% arrange(as.numeric(Position))
  teodec <- joinedteosinterep %>% filter(Chromosome == i) %>% arrange(desc(as.numeric(Position)))
  
  write_tsv(teoinc, file.path("02_Output/./teosinte_files/", paste("Teosinte_chr", i, "_increasing.txt", sep = "_")))
  write_tsv(teodec, file.path("02_Output/./teosinte_files/", paste("Teosinte_chr", i, "_decreasing.txt", sep = "_")))
}


sapply(1:10, separate_chr_files)
## Warnings NA introduced by coercion 


## Part II Visualization

#SNPs per chromosome for teosinte and maize 
joinedmaize$Chromosome <- factor(joinedmaize$Chromosome, levels = 1:10)
joinedmaizef <- joinedmaize %>%
  filter(!Chromosome %in% c("multiple", "unknown") & !is.na(Chromosome))

(plotmai <- ggplot(data = joinedmaizef) + 
  geom_bar(mapping = aes(x = Chromosome), fill = "darkgreen") +
  labs(x = "Chromosome", y = "SNP Count (Maize)", title = "Maize") +
  theme_minimal())

joinedteosinte$Chromosome <- factor(joinedteosinte$Chromosome, levels = 1:10)
joinedteosintef <- joinedteosinte %>%
  filter(!Chromosome %in% c("multiple", "unknown") & !is.na(Chromosome))

(plotteo <-ggplot(data = joinedteosintef) + 
  geom_bar(mapping = aes(x = Chromosome), fill = "yellow") +
  labs(x = "Chromosome", y = "SNP Count (Teosinte)", title = "Teosinte") +
  theme_minimal())

grid.arrange(plotmai, plotteo, ncol = 2)



(densitymai <- ggplot(joinedmaizef, aes(x = Chromosome, y = as.numeric(Position))) +
    geom_point(color = "darkgreen", alpha = 0.5) +
    labs(x = "Chromosome", y = "SNP Position", title = "Maize SNP Distribution") +
    theme_minimal())

(densityteo <- ggplot(joinedteosintef, aes(x = Chromosome, y = as.numeric(Position))) +
  geom_point(color = "yellow", alpha = 0.5) +
  labs(x = "Chromosome", y = "SNP Position", title = "Teosinte SNP Distribution") +
  theme_minimal())

grid.arrange(densitymai, densityteo, ncol = 2)


fangnew <- fang %>% select(-JG_OTU) %>% 
  pivot_longer( -Sample_ID:-Group, names_to = "SNP_ID", values_to = "Sequence")

fangnew <- fangnew %>% 
  mutate(new_sequence = ifelse(Sequence %in% c("A/A","T/T","C/C","G/G"), "Homozygous", 
                               ifelse(Sequence == "?/?", "Missing","Heterozygous")))


ggplot(fangnew, aes(x = Sample_ID, fill = new_sequence)) + geom_bar(position = "fill") + 
  theme_bw() + labs(x = "Sample ID", y = "Proportion")

ggplot(fangnew, aes(x = Group , fill = new_sequence)) + geom_bar(position = "fill") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90))+ labs(y = "Proportion")



