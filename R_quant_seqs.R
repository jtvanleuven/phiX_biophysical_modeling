###load sanger sequencing *.abif files for NNN mutagenesis and quantify sequences

library(sangerseqR)
library(stringr)
library(Biostrings)
library(seqinr)
library(reshape2)
library(kableExtra)

ref <- readDNAStringSet("rawdata/G_wt.fasta")
refaa <- Biostrings::translate(ref)
PHIX_GENETIC_CODE <- GENETIC_CODE
attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"

#function for finding codons
#loops though a list (dataframe "table") of mutated sites
#table - dataframe of sites
#fns - list of abi trace files
#ref - WT gene sequence
#minreadlength - drop read if too short
#minpeakratio - major vs. secondary peak ratio to call heterozygous
#maxmismatch - #number of allowable mismatches in the 100bp flanking the modified site
#s - start position of codon
#e - end position of codon
get_cod <- function(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e){
  #need to fix Biostrings::translate GENETIC_CODE, otherwise TTG and CTG get translated to methionine
  PHIX_GENETIC_CODE <- GENETIC_CODE
  attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"
  #set up output directories
  out.dir <- paste(dirname(fns[1]),"/failed_seqs",sep="")
  if (!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  ab1s <- list()
  reads <- list()
  hetcalls <- list()
  for (i in 1:nrow(table)){
    print(i)
    ab1s[[i]] <- readsangerseq(fns[i])
    if(length(ab1s[[i]]@primarySeq) > minreadlength) {    ####don't analyze the sequence if too short
      hetcalls[[i]] <- makeBaseCalls(ab1s[[i]], ratio = minpeakratio)
      reads[[i]] <- reverseComplement(primarySeq(ab1s[[i]]))
      #chromatogram(hetcalls[[i]], showhets=T, showcalls = "both")
      pa <- pairwiseAlignment(ref,reverseComplement(primarySeq(hetcalls[[i]])),type="global-local")
      #writePairwiseAlignments(pa)
      missum <- mismatchSummary(pa)
      if(sum(missum$pattern$position[c(max(0,(s-flankmin)):(s-1),(e+1):min((e+flankmin),max(missum$pattern$position$Position))),]$Count) < maxmismatch) {   ####don't analyze seqs that don't match reference
        qstart <- start(subject(pa))  #subject start alignment ##these are revcomp
        qend <- end(subject(pa)) #subject end alignment ##these are revcomp
        rstart <- start(pattern(pa))
        rend <- end(pattern(pa))
        chromatogram(hetcalls[[i]], trim5=length(reads[[i]])-qend, trim3=qstart, showhets=T, showcalls = "both")  ##shows region that aligns to G
        pa.het <- pairwiseAlignment(ref,reverseComplement(secondarySeq(hetcalls[[i]])), type="global-local")
        #writePairwiseAlignments(pa.het)
        missum.het <- mismatchSummary(pa.het)
        if(sum(missum.het$pattern$position[c(max(0,(s-flankmin)):(s-1),(e+1):min((e+flankmin),max(missum$pattern$position$Position))),]$Count) < maxmismatch) {   ####don't analyze seqs that don't match reference. even the alternative allele sequence should match
          get <- (s-rstart)+qstart
          indels <- sum(str_count(str_split(pa@subject,"",simplify = T)[1,1:s],"-")) ###indels mess up getting right codon. count how many and fix
          refindels <- sum(str_count(str_split(pa@pattern,"",simplify = T)[1,1:e],"-")) ###indels mess up getting right codon. count how many and fix
          primarycod <- subseq(reverseComplement(primarySeq(hetcalls[[i]])),get-indels+refindels,get+2-indels+refindels)
          secondarycod <- subseq(reverseComplement(secondarySeq(hetcalls[[i]])),get-indels+refindels,get+2-indels+refindels)
          tmp <- indel(pa)
          if(any(! s:e %in% c(tmp@insertion@unlistData@start,tmp@deletion@unlistData@start))){     #####don't analyze seqs that have indel at mut codon
            if(primarycod == secondarycod){           #####don't analyze seqs that are hetero at mut codon
              table[i,]$codon <- as.character(primarycod)
              table[i,]$aa <- as.character(Biostrings::translate(primarycod, PHIX_GENETIC_CODE))
            }else{
              table[i,]$codon <- "het"
              print("het")
            }
          }else{
            print("Failed-indel at codon")
            out.file <- paste(out.dir,paste("indel_", table$sample.names[i],".pdf",sep=""),sep="/")
            chromatogram(makeBaseCalls(ab1s[[i]], ratio = minpeakratio),filename=out.file)
          }
        }else{
          print("Failed-seq won't align")
          out.file <- paste(out.dir,paste("noalign_", table$sample.names[i],".pdf",sep=""),sep="/")
          chromatogram(makeBaseCalls(ab1s[[i]], ratio = minpeakratio),filename=out.file)
        }
      }else{
        print("Failed-seq won't align")
        out.file <- paste(out.dir,paste("noalign_", table$sample.names[i],".pdf",sep=""),sep="/")
        chromatogram(makeBaseCalls(ab1s[[i]], ratio = minpeakratio),filename=out.file)
        }
    }else{
      print("Failed-seq too short")
      out.file <- paste(out.dir,paste("tooshort_", table$sample.names[i],".pdf",sep=""),sep="/")
      chromatogram(makeBaseCalls(ab1s[[i]], ratio = minpeakratio),filename=out.file)
      }
  }
  return(table)
}

####analyze site 80 isolates 1-146
####site 80 = nuc 240 in G (2395+240=2635)
####used 2953 primer. Want ~318bp into the sequence
trim5 <- 100   ##take 100bp off 5' end
keep <- 400    ##keep 400bp after 5' trim
s <- 238    #first nucleotide of codon
e <- 240    #last nucleotide of codon
minreadlength <- 300
minpeakratio <- 0.25
maxmismatch <- 3
flankmin <- 15

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G80", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
sample.names <- str_replace(sample.names,"_80-","_G80-")  ##little fix
#####might need to change these lines depending on how the samples are named
run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,3]
sample.names <- str_split(sample.names,"_",simplify = T)[,4]
table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#a few that have heterozygous peaks - 46,80 (i<-11),88 (i<-6),90(i<-69)
#write.csv(table,file="results/g80_1-93.csv")
#write.csv(table,file="results/g80_1-147.csv")


####analyze site 8
####site 8 = nuc 24 in G (2395+24=2419)
####used 2953 primer. Want ~534bp into the sequence
setwd(wd)
setwd("G8")
trim5 <- 100   ##take 100bp off 5' end
keep <- 400    ##keep 400bp after 5' trim
s <- 22    #first nucleotide of codon
e <- 24    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches in the 200bp flanking the modified site. modified site might have 3

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G8", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
#damnit. at least 2 types of naming schemes here. don't want to rename files
sample.names <- str_replace(sample.names,"NNN-","") #helps a bit
sample.names <- str_replace(sample.names,"_NNN_","-") #helps a bit
plate.pos.1 <- str_split(sample.names,"_",simplify = T)[1:184,3]
run.date.1 <- str_split(sample.names,"_",simplify = T)[1:184,1]
sample.names.1 <- str_split(sample.names,"_",simplify = T)[1:184,4]

sample.names.2 <- str_split(sample.names,"_",simplify = T)[185:223,1]
plate.pos.2 <- str_split(sample.names,"_",simplify = T)[185:223,6]
run.date.2 <- str_split(sample.names,"_",simplify = T)[185:223,5]

plate.pos <- c(plate.pos.1,plate.pos.2)
sample.names <- c(sample.names.1,sample.names.2)
run.date <- c(run.date.1,run.date.2)

table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
#write.csv(table,file="results/g8_1-135.csv")
#write.csv(table,file="results/g8_1-223.csv")



####analyze site 74
####site 74 = nuc 222 in G (2395+222=2617)
####used 2953 primer. Want ~336bp into the sequence
trim5 <- 100   ##take 100bp off 5' end
keep <- 400    ##keep 400bp after 5' trim
s <- 220    #first nucleotide of codon
e <- 222    #last nucleotide of codon
minreadlength <- 400
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches in the 200bp flanking the modified site. modified site might have 3

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G74", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
#damnit. at least 2 types of naming schemes here. don't want to rename files
sample.names <- str_replace(sample.names,"_74","_G74") #that will fix a little
sample.names <- str_replace(sample.names,"_NNN_","-") #helps a bit more
plate.pos.1 <- str_split(sample.names,"_",simplify = T)[1:144,3]
sample.names.1 <- str_split(sample.names,"_",simplify = T)[1:144,4]
plate.pos.2 <- str_split(sample.names,"_",simplify = T)[145:193,2]
sample.names.2 <- str_split(sample.names,"_",simplify = T)[145:193,3]
sample.names <- c(sample.names.1,sample.names.2)
plate.pos <- c(plate.pos.1,plate.pos.2)
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
#write.csv(table,file="results/g74_1-199.csv")


####analyze site 15
####site 15 = nuc 45 in G (2395+45=2440)
####used 2953 primer. Want ~513bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 400    ##keep 400bp after 5' trim
s <- 43    #first nucleotide of codon
e <- 45    #last nucleotide of codon
minreadlength <- 400
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches in the 200bp flanking the modified site. modified site might have 3

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G15", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
#damnit. at least 2 types of naming schemes here. don't want to rename files
sample.names <- str_replace(sample.names,"_15-","_G15-") #that will fix a little
sample.names <- str_replace(sample.names,"_NNN_","-") #helps a bit more
plate.pos <- str_split(sample.names,"_",simplify = T)[,2]
sample.names <- str_extract(sample.names,"G15-[0-9a-z]+")
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
#write.csv(table,file="results/g15_1-200.csv")


####analyze site 15 NNK
####site 15 = nuc 45 in G (2395+45=2440)
####used 2953 primer. Want ~513bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 400    ##keep 400bp after 5' trim
s <- 43    #first nucleotide of codon
e <- 45    #last nucleotide of codon
minreadlength <- 400
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking the modified site. modified site might have 3

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G15_NNK", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
#damnit. at least 2 types of naming schemes here. don't want to rename files
sample.names <- str_replace(sample.names,"NNK-","") #that will fix a little
sample.names <- str_replace(sample.names,"_NNK-","-") #helps a bit more
sample.names <- str_replace(sample.names,"G15_","G15-") #helps a bit more
plate.pos.1 <- str_split(sample.names,"_",simplify = T)[1:50,3]
run.date.1 <- str_split(sample.names,"_",simplify = T)[1:50,1]
sample.names.1 <- str_split(sample.names,"_",simplify = T)[1:50,4]

sample.names.2 <- str_split(sample.names,"_",simplify = T)[51:100,1]
plate.pos.2 <- str_split(sample.names,"_",simplify = T)[51:100,6]
run.date.2 <- str_split(sample.names,"_",simplify = T)[51:100,5]

plate.pos <- c(plate.pos.1,plate.pos.2)
sample.names <- c(sample.names.1,sample.names.2)
run.date <- c(run.date.1,run.date.2)

table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
#write.csv(table,file="results/g15_NNK_1-100.csv")



####analyze site 123 isolates 1-93
####analyze site 123 isolates 1-144
####site 123 = nuc 369 in G (2395+369=2764)
####used 2953 primer. Want ~189bp into the sequence
trim5 <- 50   ##take 50bp off 5' end
keep <- 500    ##keep 500bp after 5' trim
s <- 367    #first nucleotide of codon
e <- 369    #last nucleotide of codon
minreadlength <- 300
minpeakratio <- 0.3  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site
flankmin <- 15

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G123", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
View(sample.names)
#####might need to change these lines depending on how the samples are named
#damnit. at least 2 types of naming schemes here. don't want to rename files
plate.pos.1 <- str_split(sample.names,"_",simplify = T)[1:48,2]
run.date.1 <- str_split(sample.names,"_",simplify = T)[1:48,1]
sample.names.1 <- str_split(sample.names,"_",simplify = T)[1:48,3]

sample.names.2 <- str_split(sample.names,"_",simplify = T)[49:186,1]
plate.pos.2 <- str_split(sample.names,"_",simplify = T)[49:186,6]
run.date.2 <- str_split(sample.names,"_",simplify = T)[49:186,5]

plate.pos <- c(plate.pos.1,plate.pos.2)
sample.names <- c(sample.names.1,sample.names.2)
run.date <- c(run.date.1,run.date.2)

table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#a few that have heterozygous peaks - 46,80 (i<-11),88 (i<-6),90(i<-69)
#write.csv(table,file="results/g123_1-94.csv")
#write.csv(table,file="results/g123_1-144.csv")
#write.csv(table,file="results/g123_1-186.csv")
write.csv(table,file="results/g123.csv")



####analyze site 8_NNK isolates 1-190 (picked by stinger)
####site 8 = nuc 24 in G (2395+24=2419)
####used 2953 primer. Want ~534bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 22    #first nucleotide of codon
e <- 24    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G8_NNK", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"G8_NNK","G8-NNK")
run.date <- str_split(sample.names,"_",simplify = T)[,5]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
table <- table[!table$codon=="TTC",]
#write.csv(table,file="results/g8_nnk_1-190.csv")


####analyze site 106 isolates 1-90
####reanalyze site 106 isolates 1-192
####site 106 = nuc 318 in G (2395+318=2713)
####used 2953 primer. 2953-2713 Want ~240bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 316    #first nucleotide of codon
e <- 318    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G106", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"G106_NNN","G106")
run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,3]
sample.names <- str_split(sample.names,"_",simplify = T)[,4]
table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
###table <- table[!table$codon=="TTC",]   ###easy way to remove one bad one
#write.csv(table,file="results/g106_1-95.csv")
#write.csv(table,file="results/g106_1-192.csv")



####analyze site 117 isolates 1-96
####analyze site 117 isolates 1-192
####site 117 = nuc 351 in G (2395+351=2746)
####used 2953 primer. 2953-2746 Want ~207bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 349    #first nucleotide of codon
e <- 351    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G117", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"G117_NNN","G117")
run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,3]
sample.names <- str_split(sample.names,"_",simplify = T)[,4]
table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
###table <- table[!table$codon=="TTC",]   ###easy way to remove one bad one
#write.csv(table,file="results/g117_1-96.csv")
#write.csv(table,file="results/g117_1-192.csv")



####analyze site 43 isolates 1-96
####analyze site 43 isolates 1-192
####site 43 = nuc 129 in G (2395+129=2524)
####used 2953 primer. 2953-2524 Want ~429bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 127    #first nucleotide of codon
e <- 129    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G43", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"G43_NNN","G43")
run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,3]
sample.names <- str_split(sample.names,"_",simplify = T)[,4]
table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
###table <- table[!table$codon=="TTC",]   ###easy way to remove one bad one
#write.csv(table,file="results/g43_1-96.csv")
#write.csv(table,file="results/g43_1-192.csv")


####analyze site 41 isolates 1-96
####analyze site 41 isolates 1-144
####site 41 = nuc 123 in G (2395+123=2518)
####used 2953 primer. 2953-2518 Want ~435bp into the sequence
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 121    #first nucleotide of codon
e <- 123    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site
flankmin <- 50

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G41", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"G41_NNN","G41")
sample.names <- str_replace(sample.names,"_\\d\\d\\d_","_")
run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,2]
sample.names <- str_split(sample.names,"_",simplify = T)[,3]
table <- cbind(run.date,plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
#table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
###table <- table[!table$codon=="TTC",]   ###easy way to remove one bad one
#write.csv(table,file="results/g41_1-96.csv")
#write.csv(table,file="results/g41_1-144.csv")


####analyze site 88 isolates 1-
####site 88 = nuc 264 in G (2395+264=2659)
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 262    #first nucleotide of codon
e <- 264    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site
flankmin <- 50

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G88", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX2953","phiX_2953")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,7]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#write.csv(table,file="results/g88_1-190.csv")
#write.csv(table,file="results/g88_1-234.csv")
write.csv(table,file="results/g88.csv")


####analyze site 119 isolates 1-94
####analyze site 119 isolates 1-190
####site 119*3 = nuc 357 in G (2395+357=2752)
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 355    #first nucleotide of codon
e <- 357    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site
flankmin <- 50

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G119", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#write.csv(table,file="results/g119_1-94.csv")
#write.csv(table,file="results/g119_1-190.csv")
#write.csv(table,file="results/g119_1-235.csv")
write.csv(table,file="results/g119.csv")

####analyze site 125 isolates 1-190
####site 125*3 = nuc 375 in G (2395+375=2770)
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 373    #first nucleotide of codon
e <- 375    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G125", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
#####might need to change these lines depending on how the samples are named
#sample.names <- str_replace(sample.names,"G41_NNN","G41")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,7]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, s, e)
#write.csv(table,file="results/g125_1-190.csv")

####analyze site 129 isolates 1-144
####site 129*3 = nuc 387 in G (2395+387=2782)
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 385    #first nucleotide of codon
e <- 387    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches in the 200bp flanking(50 bp either side) the modified site
flankmin <- 20

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G129", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
#sample.names <- str_replace(sample.names,"G41_NNN","G41")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,7]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#write.csv(table,file="results/g129_1-94.csv")
#write.csv(table,file="results/g129_1-144.csv")


####analyze site 145 isolates 1-127
####site 148*3 = nuc 444 in G (2395+444=2839)
trim5 <- 10   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 442    #first nucleotide of codon
e <- 444    #last nucleotide of codon
minreadlength <- 500
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 10 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G148", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
#sample.names <- str_replace(sample.names,"G41_NNN","G41")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,7]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#write.csv(table,file="results/g148_1-94.csv")
#write.csv(table,file="results/g148_1-127.csv")



####analyze site 4 isolates 1-127
####site 4*3 = nuc 12 in G (2395+12=2407)
trim5 <- 10   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 10    #first nucleotide of codon
e <- 12    #last nucleotide of codon
minreadlength <- 600
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 10 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G4", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
#run.date <- str_split(sample.names,"_",simplify = T)[,5]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
#write.csv(table,file="results/g4_1-127.csv")
#write.csv(table,file="results/g4_1-173.csv")
write.csv(table,file="results/g4.csv")


####analyze site 72 isolates 1-159
####site 72*3 = nuc 216 in G (2395+216=2611)
trim5 <- 10   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 214    #first nucleotide of codon
e <- 216    #last nucleotide of codon
minreadlength <- 500
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 20 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G72", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g72_1-159.csv")


###analyze site 3 isolates 1-X
####site 3*3 = nuc 9 in G (2395+9=2404)
trim5 <- 10   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 7    #first nucleotide of codon
e <- 9    #last nucleotide of codon
minreadlength <- 600
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 5 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G3", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
sample.names <- str_replace(sample.names, "_NNN_", "-")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g3.csv")


####analyze site 11 isolates 1-X
####site 11*3 = nuc 33 in G (2395+33=2428)
trim5 <- 10   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 31    #first nucleotide of codon
e <- 33    #last nucleotide of codon
minreadlength <- 600
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 5 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G11", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
sample.names <- str_replace(sample.names, "_NNN-", "-")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g11.csv")

####analyze site 145 isolates 1-X
####site 45*3 = nuc 135 in G (2395+135=2530)
trim5 <- 10   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 133    #first nucleotide of codon
e <- 135    #last nucleotide of codon
minreadlength <- 600
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 5 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G45", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
sample.names <- str_replace(sample.names, "_NNN-", "-")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g45.csv")


####analyze site 128 isolates 1-X
####site 128*3 = nuc 384 in G (2395+384=2779)
trim5 <- 5   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 382    #first nucleotide of codon
e <- 384    #last nucleotide of codon
minreadlength <- 600
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 5 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G128", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
sample.names <- str_replace(sample.names, "_NNN-", "-")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g128.csv")



###analyze site 4 isolates 1-X
####site 3*4 = nuc 12 in G (2395+12=2407)
trim5 <- 10   ##take 100bp off 5' end
keep <- 700    ##keep 400bp after 5' trim
s <- 10    #first nucleotide of codon
e <- 12    #last nucleotide of codon
minreadlength <- 600
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 4 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 5 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G4_new_NNN_primers", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
sample.names <- str_replace(sample.names, "_NNN-", "-")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g4_reorder.csv")


###analyze site 125 isolates 1-X
####site 3*125 = nuc 375 in G (2395+375=2770)
trim5 <- 50   ##take 100bp off 5' end
keep <- 650    ##keep 400bp after 5' trim
s <- 373    #first nucleotide of codon
e <- 375    #last nucleotide of codon
minreadlength <- 550
minpeakratio <- 0.25  #minimum minor peak to major peak ratio to call heterozygous site
maxmismatch <- 3 #number of allowable mismatches with "flankmin" of the modified site
flankmin <- 5 #number of bases 5' and 3' of the modified codon needed to match btween read and seq 

wt.cod <- subseq(ref,s,e)
fns <- sort(list.files("rawdata/sanger_reads/G125_new_NNN_primers", full.names = TRUE))
fns <- fns[grepl("ab1", fns)]
sample.names <- sapply(strsplit(fns, ".ab1"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])
head(sample.names)
#####might need to change these lines depending on how the samples are named
sample.names <- str_replace(sample.names,"phiX_2953R","phiX2953R")
sample.names <- str_replace(sample.names, "_NNN-", "-")
#run.date <- str_split(sample.names,"_",simplify = T)[,1]
plate.pos <- str_split(sample.names,"_",simplify = T)[,6]
sample.names <- str_split(sample.names,"_",simplify = T)[,1]
table <- cbind(plate.pos,sample.names)
table <- as.data.frame(table)
table$codon <- NA
table$aa <- NA
table <- get_cod(table, fns, ref, minreadlength, minpeakratio, maxmismatch, flankmin, s, e)
sum(table(table$aa))/nrow(table)
#write.csv(table,file="results/g125_reorder.csv")




#####a little random analysis
#wts G8-R , G15-S , G74-S, G80-A, G123-V, G106-A, G117-D, G43-Y
g3 <- read.csv("rawdata/sanger_reads/g3.csv", header=T, row.names = 1)
g4 <- read.csv("rawdata/sanger_reads/g4.csv", header=T, row.names = 1)
g4_re <- read.csv("rawdata/sanger_reads/g4new.csv", header=T, row.names = 1)
g8 <- read.csv("rawdata/sanger_reads/g8.csv", header=T, row.names = 1)
g8_NNK <- read.csv("rawdata/sanger_reads/g8_nnk.csv", header=T, row.names = 1)
g11 <- read.csv("rawdata/sanger_reads/g11.csv", header=T, row.names = 1)
g15 <- read.csv("rawdata/sanger_reads/g15.csv", header=T, row.names = 1)
g15_NNK <- read.csv("rawdata/sanger_reads/g15_NNK.csv", header=T, row.names = 1)
g41 <- read.csv("rawdata/sanger_reads/g41.csv", header=T, row.names = 1)
g43 <- read.csv("rawdata/sanger_reads/g43.csv", header=T, row.names = 1)
g45 <- read.csv("rawdata/sanger_reads/g45.csv", header=T, row.names = 1)
g72 <- read.csv("rawdata/sanger_reads/g72.csv", header=T, row.names = 1)
g74 <- read.csv("rawdata/sanger_reads/g74.csv", header=T, row.names = 1)
g80 <- read.csv("rawdata/sanger_reads/g80.csv", header=T, row.names = 1)
g81 <- read.csv("rawdata/sanger_reads/g81.csv", header=T, row.names = 1)
g88 <- read.csv("rawdata/sanger_reads/g88.csv", header=T, row.names = 1)
g106 <- read.csv("rawdata/sanger_reads/g106.csv", header=T, row.names = 1)
g117 <- read.csv("rawdata/sanger_reads/g117.csv", header=T, row.names = 1)
g119 <- read.csv("rawdata/sanger_reads/g119.csv", header=T, row.names = 1)
g123 <- read.csv("rawdata/sanger_reads/g123.csv", header=T, row.names = 1)
g125 <- read.csv("rawdata/sanger_reads/g125.csv", header=T, row.names = 1)
g125_re <- read.csv("rawdata/sanger_reads/g125new.csv", header=T, row.names = 1)
g128 <- read.csv("rawdata/sanger_reads/g128.csv", header=T, row.names = 1)
g129 <- read.csv("rawdata/sanger_reads/g129.csv", header=T, row.names = 1)
g148 <- read.csv("rawdata/sanger_reads/g148.csv", header=T, row.names = 1)



#how many codons?
all_cods <- rbind(g3[,names(g3) %in% c("sample.names", "codon")],
                  g4[,names(g4) %in% c("sample.names", "codon")],
                  g4_re[,names(g4_re) %in% c("sample.names", "codon")],
                  g8[,names(g8) %in% c("sample.names", "codon")],
                  g11[,names(g11) %in% c("sample.names", "codon")],
                  g15[,names(g15) %in% c("sample.names", "codon")],
                  g41[,names(g41) %in% c("sample.names", "codon")],
                  g43[,names(g43) %in% c("sample.names", "codon")],
                  g45[,names(g45) %in% c("sample.names", "codon")],
                  g72[,names(g72) %in% c("sample.names", "codon")],
                  g74[,names(g74) %in% c("sample.names", "codon")],
                  g80[,names(g80) %in% c("sample.names", "codon")],
                  g81[,names(g81) %in% c("sample.names", "codon")],
                  g88[,names(g88) %in% c("sample.names", "codon")],
                  g106[,names(g106) %in% c("sample.names", "codon")],
                  g117[,names(g117) %in% c("sample.names", "codon")],
                  g119[,names(g119) %in% c("sample.names", "codon")],
                  g123[,names(g123) %in% c("sample.names", "codon")],
                  g125[,names(g125) %in% c("sample.names", "codon")],
                  g125_re[,names(g125_re) %in% c("sample.names", "codon")],
                  g128[,names(g128) %in% c("sample.names", "codon")],
                  g129[,names(g129) %in% c("sample.names", "codon")],
                  g148[,names(g148) %in% c("sample.names", "codon")])

tot <- nrow(all_cods)
print(tot)
failed <- unname(table(is.na(all_cods$codon))[2])
print(failed)
hets <- unname(table(all_cods$codon == "het")[2])
print(hets)
cods <- length(table(all_cods$codon))
all_cods_old <- all_cods
ktable <- data.frame(c(tot,failed,hets,cods))
names(ktable) <- "Seq count"
row.names(ktable) <- c("Total", "Failed", "Hets", "Unique codons")

kable(ktable, "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), position="center", font_size = 16, options(knitr.kable.NA = '-')) %>%
  column_spec(1, bold = T) %>%
  column_spec(2, width = "30em")

all_cods <- na.omit(all_cods[!all_cods$codon == "het",])
all_cods$isolate <- paste(str_split(all_cods$sample.names,"-",simplify = T)[,1],all_cods$codon,sep="_")
length(unique(all_cods$isolate)) #585
#write.csv(all_cods, "results/all_cods.csv", quote = F, row.names = F)  ##has both redo and original genotyping results

###get me some basic stats
filt <- data.frame(table(str_split(all_cods$sample.names, "-", simplify = T)[,1]))
unfilt <- data.frame(table(str_split(all_cods_old$sample.names, "-", simplify = T)[,1]))
ktable <- merge(filt,unfilt,by="Var1")
ktable$prop <- ktable$Freq.x/ktable$Freq.y
ktable$prop <- round(ktable$prop, digits = 2)
ktable.t <- data.frame(t(ktable), stringsAsFactors = FALSE)
colnames(ktable.t) <- t(ktable$Var1)
ktable.t <- ktable.t[-1,]
row.names(ktable.t) <- c("Good", "All", "Freq")


kable(ktable.t, "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), position="center", font_size = 16, options(knitr.kable.NA = '-')) %>%
  column_spec(1, bold = T) %>%
  column_spec(2:17, width = "30em")

ktable2 <- ktable
names(ktable2) <- c("Var1","NumReads","RawReads","PassFreq")
ktable2$Var1 <- as.numeric(str_replace(ktable2$Var1,"G",""))
lookup <- str_split(as.character(refaa),"")
ktable2$ref <- lookup[[1]][ktable2$Var1]
ktable2 <- ktable2[order(ktable2$Var1),]
ktable2$Residue <- paste(ktable2$ref,ktable2$Var1,sep="")
ktable2$NumCods <- NA
ktable2$NumAAs <- NA
all.cat <- all
all.cat$site <- str_replace("re","")
for(i in 1:nrow(ktable2)){
  get <- paste("G",ktable2[i,]$Var1, sep="")
  cod.s <- all_cods[str_split(all_cods$sample.names,"-",simplify = T)[,1]==get,]
  ktable2[i,]$NumCods <- length(unique(cod.s$isolate))
  aa.s <- all.cat[all.cat$site == get,]
  ktable2[i,]$NumAAs <- length(unique(aa.s$aa))
}

kable(ktable2[,c("Residue","NumReads","NumCods","NumAAs")], "html", row.names = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), position="center", font_size = 12, options(knitr.kable.NA = '-')) %>%
  column_spec(1, bold = T) %>%
  column_spec(2:4, width = "30em")



aa_g3 <- cbind(data.frame(table(g3$aa)),rep("G3",length(table(g3$aa))))
names(aa_g3) <- c("aa","freq","site")
aa_g4 <- cbind(data.frame(table(g4$aa)),rep("G4",length(table(g4$aa))))
names(aa_g4) <- c("aa","freq","site")
aa_g4_re <- cbind(data.frame(table(g4_re$aa)),rep("G4re",length(table(g4_re$aa))))
names(aa_g4_re) <- c("aa","freq","site")
aa_g8 <- cbind(data.frame(table(g8$aa)),rep("G8",length(table(g8$aa))))
names(aa_g8) <- c("aa","freq","site")
aa_g11 <- cbind(data.frame(table(g11$aa)),rep("G11",length(table(g11$aa))))
names(aa_g11) <- c("aa","freq","site")
aa_g15 <- cbind(data.frame(table(g15$aa)),rep("G15",length(table(g15$aa))))
names(aa_g15) <- c("aa","freq","site")
#aa_g15_NNK <- cbind(data.frame(table(g15_NNK$aa)),rep("G15_NNK",length(table(g15_NNK$aa))))
#names(aa_g15_NNK) <- c("aa","freq","site")
aa_g41 <- cbind(data.frame(table(g41$aa)),rep("G41",length(table(g41$aa))))
names(aa_g41) <- c("aa","freq","site")
aa_g43 <- cbind(data.frame(table(g43$aa)),rep("G43",length(table(g43$aa))))
names(aa_g43) <- c("aa","freq","site")
aa_g45 <- cbind(data.frame(table(g45$aa)),rep("G45",length(table(g45$aa))))
names(aa_g45) <- c("aa","freq","site")
aa_g72 <- cbind(data.frame(table(g72$aa)),rep("G72",length(table(g72$aa))))
names(aa_g72) <- c("aa","freq","site")
aa_g74 <- cbind(data.frame(table(g74$aa)),rep("G74",length(table(g74$aa))))
names(aa_g74) <- c("aa","freq","site")
aa_g80 <- cbind(data.frame(table(g80$aa)),rep("G80",length(table(g80$aa))))
names(aa_g80) <- c("aa","freq","site")
aa_g81 <- cbind(data.frame(table(g81$aa)),rep("G81",length(table(g81$aa))))
names(aa_g81) <- c("aa","freq","site")
aa_g88 <- cbind(data.frame(table(g88$aa)),rep("G88",length(table(g88$aa))))
names(aa_g88) <- c("aa","freq","site")
aa_g106 <- cbind(data.frame(table(g106$aa)),rep("G106",length(table(g106$aa))))
names(aa_g106) <- c("aa","freq","site")
aa_g117 <- cbind(data.frame(table(g117$aa)),rep("G117",length(table(g117$aa))))
names(aa_g117) <- c("aa","freq","site")
aa_g119 <- cbind(data.frame(table(g119$aa)),rep("G119",length(table(g119$aa))))
names(aa_g119) <- c("aa","freq","site")
aa_g123 <- cbind(data.frame(table(g123$aa)),rep("G123",length(table(g123$aa))))
names(aa_g123) <- c("aa","freq","site")
aa_g125 <- cbind(data.frame(table(g125$aa)),rep("G125",length(table(g125$aa))))
names(aa_g125) <- c("aa","freq","site")
aa_g125_re <- cbind(data.frame(table(g125_re$aa)),rep("G125re",length(table(g125_re$aa))))
names(aa_g125_re) <- c("aa","freq","site")
aa_g128 <- cbind(data.frame(table(g128$aa)),rep("G128",length(table(g128$aa))))
names(aa_g128) <- c("aa","freq","site")
aa_g129 <- cbind(data.frame(table(g129$aa)),rep("G129",length(table(g129$aa))))
names(aa_g129) <- c("aa","freq","site")
aa_g148 <- cbind(data.frame(table(g148$aa)),rep("G148",length(table(g148$aa))))
names(aa_g148) <- c("aa","freq","site")


all <- rbind(aa_g3,aa_g4,aa_g4_re,aa_g8,aa_g11,aa_g15,aa_g41,aa_g43,aa_g45,aa_g72,aa_g74,aa_g80,aa_g81,aa_g88,aa_g106,aa_g117,aa_g119,aa_g123,aa_g125,aa_g125_re,aa_g128,aa_g129,aa_g148)
list <- names(table(str_split(all_cods_old$sample.names, "-", simplify = T)[,1]))
list[! list %in% names(table(all$site))]

all.wide <- data.frame(t(dcast(all, aa ~ site, value.var="freq")))
names(all.wide) <- as.character(as.matrix(all.wide[1,]))
all.wide <- all.wide[-1,]

PHIX_GENETIC_CODE <- GENETIC_CODE
attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"
aas <- unique(as.character(PHIX_GENETIC_CODE)) 
aas.miss <- aas[! aas %in% names(all.wide)]
miss <- data.frame(matrix(nrow=nrow(all.wide), ncol=length(aas.miss)))
names(miss) <- aas.miss

all.wide <- cbind(all.wide,miss)
#write.csv(all.wide, "results/protein_G_mutagenesis_AA_counts.csv")
