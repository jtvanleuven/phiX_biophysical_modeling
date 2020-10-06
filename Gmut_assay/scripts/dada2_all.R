
#######################################################
# G Mutants Growth Assay Sequence Analysis Workflow
#######################################################

#### SUMMARY OF THE WORKFLOW ####
#### command-line:
# 1. fastp to QC
# 2. FLASh to merge forwqrd and reverse reads
# 3. samtools to extract mapped reads with CIGAR 528M (ensure identical length)
#     samtools view -h in.sam | awk '$6 == "528M"' > out.sam
# 4. Reattach header:
#     samtools view -H file.sam > sam_head.txt
#     cat sam_head.txt headerless.sam > out.sam
# 5. Convert to fastq: 
#     picard SamToFastq I=input.sam F=output.fastq
#    Make sure no Ns in fastq (dada will throw an error)
# 
#### R
# 6. run dada2 on the fastq from step 5 (R code below)
# 7. export to fasta (uniquesToFasta)
# 8. map fasta to sam: bwa mem ../../../phix_anc.fasta G1G2_t35.fasta > G1G2t35_dada.sam
# 9*. Remove header from sam files and read it into R as table again, for matching muts and counts
# 10. Downstream analysis in R (this file and post_dada2_analysis.Rmd)


# FYI: step 5 -- files are in ~/Projects/phix_host/PrelimAssay/data/processed_data/filtered/sam/

# * disregard step 9. Included skip lines option in read.table function (20200629) so no need to
#             manually remove the headers before reading into R.


#------------------------- dada2 part in this section done in server -----------------------------
#setwd("~/phix_prelim/")  
setwd("~/phix_zoe/")

library(dada2)
library(stringr)

#path <- "~/Projects/phix_host/PrelimAssay/dada2_all/"
#path <- "~/phix_prelim/"
path <- "~/phix_zoe/"
#files <- list.files(paste(path, "input_fastq/", sep = ""))
files <- list.files(paste(path, "converted_fastq/", sep = ""))
#files.in <- paste("input_fastq/", files, sep = "")
files.in <- paste("converted_fastq/", files, sep = "")
#samples <- str_split(files, ".extracted.", simplify = T)[,1]
samples <- str_split(files, "_converted.", simplify = T)[,1]

for (i in 3:length(samples)){
  err.i <- learnErrors(files.in[i], multithread=TRUE)
  dds.i <- dada(files.in[i], err=err.i, multithread=TRUE, OMEGA_A = 1e-10, OMEGA_C = 1e-10, 
                   DETECT_SINGLETONS = TRUE)
  uniquesToFasta(getUniques(dds.i), paste(path, "uniquesToFasta_out/", samples[i], ".fasta", sep = ""))
}
  
#------------------------------------------------------------------------------------------------------



# Map fasta from dada2:
#     bwa mem ../../../phix_anc.fasta G1G2_t35.fasta

#------------desregard this section--------------------------------------------
# then remove first 2 lines (header) so that sam can be read in as a table in R:
#     sed -i '1,2d' file.sam 
#     (-i deletes it from the source file. Otherwise sed sends to stout)
# For multiple files:
# for i in *.sam; do
#      sed -i '1,2d' $i
# done
# -i throwing error. do redirect (wrote a script "sed_sam.sh")
#-------------------------------------------------------------------------------


### Edited final sam in sed_sam dir

#### Separate script fn_dadasam.R contains 2 functions:
#### get_glib() to read in JT's primers file to get prot G muts
#### analyze_dadasam() to loop through sam files and get counts, muts info

library(tidyverse)
library(Biostrings)
source("scripts/fn_dadasam.R")

glib <- get_glib()

# Below code will read through sam files and organize dada2 result in a dataframe 
# with each row for each codon change

#dir <- "./data/sed_sam/"
dir <- "data/final_sam/" # with skip lines option added in fn_dadasam.R, use this directory instead
files <- list.files(dir)
samples <- str_split(files, "\\.", simplify = TRUE)[,1]
outlist <- list()
for (k in 1:length(samples)){
  outlist[[k]] <- analyze_dadasam(paste(dir, files[k], sep = ""))
}
names(outlist) <- samples
yesol_df <- data.frame()
for (i in 1:length(samples)){
  outlist[[i]]$mut_info$sample <- samples[i]
  yesol_df <- rbind(yesol_df, outlist[[i]]$mut_info)
}
yesol_df <- select(yesol_df, sample, everything())
write.csv(yesol_df, "data/dada2result_yesol.csv")

# Zoe's samples 
#dir <- "./data/sed_sam_zoe/"
dir <- "data/final_sam_20200630/"
files <- list.files(dir)
samples <- str_split(files, "\\.", simplify = TRUE)[,1]
zoelist <- list()
for (i in 1:length(samples)){
  zoelist[[i]] <- analyze_dadasam(paste(dir, files[i], sep = ""))
}
names(zoelist) <- samples
zoe_df <- data.frame()
for (i in 1:length(zoelist)){
  zoelist[[i]]$mut_info$sample <- samples[i]
  zoe_df <- rbind(zoe_df, zoelist[[i]]$mut_info)
}
zoe_df <- select(zoe_df, sample, everything())
write.csv(zoe_df, "./data/dada2result_zoe.csv")


#----------------------------------------------------------------------------

# Below code collapses mutations in a sequence into one
# i.e., each row will be a sequence var instead of a mut

yesol_df <- read.csv("data/dada2result_yesol.csv", stringsAsFactors = F)
yesol_df <- yesol_df[,-1] # First col in csv is just index numbers.
zoe_df <- read.csv("data/dada2result_zoe.csv", stringsAsFactors = F)
zoe_df <- zoe_df[,-1]

# df <- yesol_df
df <- zoe_df
df$ID <- paste(df$sample, df$var, sep = ";") #create a unique id for each var 

id_v <- unique(df$ID)
df2 <- as.data.frame(matrix(0, nrow = length(id_v), ncol = 6))
colnames(df2) <- c("sample","var","AAsub","count","nmut_tot","nmut_offtarget")
for (i in 1:length(id_v)){
  id_i <- id_v[i]
  df_i <- filter(df, ID == id_i)
  df2$sample[i] <- df_i$sample[1]
  df2$var[i] <- df_i$var[1]
  df2$count[i] <- df_i$count[1]
  df2$nmut_tot[i] <- max(df_i$nmut)
  df2$nmut_offtarget[i] <- df2$nmut_tot[i] - sum(df_i$in_target)
  for (k in 1:nrow(df_i)){
    df_i$AAsub[k] <- sub(">", df_i$site[k], df_i$AA_change[k])
  }
  df2$AAsub[i] <- paste(df_i$AAsub, collapse = ",")
}
# yesol_df2 <- df2
zoe_df2 <- df2

#write.csv(yesol_df2, "data/dada2result_yesol2.csv")
write.csv(zoe_df2, "data/dada2result_zoe2.csv")


