#Code to analyze starting codon frequencies in the NNN mutagenesis before transforming into E. coli cells


library(stringr)
library(Biostrings)
library(xlsx)
library(Rsamtools)
library(Rmisc)

PHIX_GENETIC_CODE <- GENETIC_CODE
attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"

##The initial set of Illumina sequencing was done for sites 4, 8, 15, 41, 43, 72, 74, 80, 81, 88, 106, 117, 119, 123, 125, 129, and 148.
##The ligation mixes were PCR amplified (no triplicate PCR) and sequenced by MiSeq at the beginning of August 2018 (080218)

##These were used to identify skewed initial starting pools for sites 4 and 125, and new primers were made. 
##We then sequenced the new ligation mixes and added site 3. These PCRs were done in triplicate.
##sequencing was done at the end of April 2019 (042219 or 05XX19)


count.cods <- function(aln,s,e,res){
  PHIX_GENETIC_CODE <- GENETIC_CODE
  attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"
  obs <- NULL
  for(i in 1:length(aln$qname)){
    seq <- aln$seq[i]
    qual <- aln$qual[i]
    cig <- aln$cigar[i]
    pos <- aln$pos[i]
    if(!is.na(aln$mapq[i])){    #some reads have no mapping quality
      if(aln$mapq[i] > 50){     #minimum mapping quality
        cig2 <- str_split(cig,"[A-Z]",simplify = T)
        cig3 <- str_split(cig,"\\d+",simplify = T)[,-1]
        if(cig2[1] == '11' && cig3[1] == 'I'){    ###fix problem with reads hanging over edge of reference sequence. Currently designated by "I" in cigar string
              cig <- str_replace(cig, '11I', '11S')
              cig3[1] <- 'S'
        }
        if(s > pos && s < pos+sum(as.numeric(cig2[which(cig3=="M")]))){  ##make sure the codon of interest is within the mapping range
          if(all(!str_detect(cig, c("I","D","P")))){   #ignore reads with indels
            if(cig3[1] == "S" | cig3[1] == "H"){
              s2 <- (s - pos) + as.numeric(cig2[1]) + 1
            #}else if(cig3[1] == "I"){
              #s2 <- (s - aln$pos[i]) + as.numeric(cig2[1]) + 1
            }else if(cig3[1] == "M"){
              s2 <- s - (pos - 1)
            }else{
              print(paste("unreconginzed cigar:", cig, sep=" "))
            }
            e2 <- s2+2
            cod <- substr(as.character(seq),s2,e2)
            startcod <- substr(as.character(seq),s2-s+1,s2-s+3)
            #if(startcod == 'ATG'){                                      ##only useful for forward reads
              #paste(substr(as.character(seq),s2-10,s2-1),substr(as.character(seq),s2,e2),substr(as.character(seq),e2+1,e2+10),sep=" ")
              cod.qual <- substr(as.character(qual),s2,e2)
              #paste(substr(as.character(qual),s2-10,s2-1),substr(as.character(qual),s2,e2),substr(as.character(qual),e2+1,e2+10),sep=" ")
              #as.numeric(as.character(charToRaw(cod.qual))) - 33 + 24 #not really sure why R doesn't get ascii right
              #qual[cod]
              #aln$pos[i]
              obs <- c(obs,cod)
            #}else{
            #  print(paste("Not ATG at ", i, sep = ''))
            #}
          }
        }
      }
    }
  }
  obs.tab <- table(obs)
  #cod.tab[j,] <- as.numeric(obs.tab[names(PHIX_GENETIC_CODE)])
  cod.tab <- as.numeric(obs.tab[names(PHIX_GENETIC_CODE)])
  #print(paste("done with ",tab$id[j],sep=""))
  return(obs.tab)
}


##on linux for demultiplexing reads
# module load python
# source /mnt/lfs2/grcuser/opt/dbcAmpVenv/bin/activate
# splitReadsBySample.py -1 phix_41_43_106_R1.fastq -2 phix_41_43_106_R2.fastq
# splitReadsBySample.py -1 May2019_phiX_R1.fastq.gz -2 May2019_phiX_R2.fastq.gz
# run_variant_calling.sh

ref <- readDNAStringSet("rawdata/G_wt.fasta")



####process reads from 042219
sam.tab <- read.xlsx("rawdata/SubmissionSheet_Amplicon_v3.0_042219_jtv.xls", sheetIndex = 1)
sam.tab <- sam.tab[c(16,17,18,151:156),1:14]
names(sam.tab) <- as.character(as.matrix(sam.tab[1,]))
sam.tab <- sam.tab[-1,]

fns <- sort(list.files("rawdata/illumina_reads/042219/bam/", full.names = TRUE))
fnFs <- fns[grepl("sorted.bam.bai", fns)]
sample.names <- sapply(strsplit(fnFs, ".aligned"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])

#g41 <- data.frame(sample.names[str_detect(sample.names, "g41_t")])
tab <- data.frame(sample.names)
names(tab) <- "id"
tab$id <- str_replace(tab$id,"-","_")
tab$id <- str_replace(tab$id,"_redo","-redo")
tab$time <- str_split(tab$id,"_",simplify = T)[,2]
tab$rep <- str_split(tab$id,"_",simplify = T)[,3]
tab$site <- str_split(tab$id,"_",simplify = T)[,1]
tab$site <- str_extract(tab$site,"\\d+")
tab$id <- str_replace(tab$id,"_R1","")
tab$id <- str_replace(tab$id,"-redo","_redo")
#tab[which(str_detect(tab$id, "NTC")),]$time <- 90 
#tab[which(str_detect(tab$id, "NTC")),]$rep <- c(2,3,1)


##could not find a good tool to measure codon frequencies over time.
#have to write myself
#####try with bam files
###code works with single mix (not multiple sites mixed together)
cod.tab <- data.frame(matrix(nrow=nrow(tab),ncol=64))
names(cod.tab) <- names(PHIX_GENETIC_CODE)
slow <- F
if(slow){
  for(j in 1:nrow(tab)){
    tab$id[j]
    bam.sorted <- BamFile(paste("rawdata/illumina_reads/042219/bam/",tab$id[j],".aligned.sorted.bam",sep=""))
    #seqinfo(bam.sorted)
    #names(bam.sorted)
    #quickBamFlagSummary(bam.sorted)
    aln <- scanBam(bam.sorted)
    aln <- aln[[1]]
    #names(aln)
    #lapply(aln, function(xx) xx[1])
    cigs <- table(aln$cigar)
    cigs[cigs > length(aln$qname)*0.001]
    #filter alinged reads by
    #1 mapq
    #2 isize
    res <- as.numeric(tab$site[j])
    s <- res*3-2
    e <- res*3
    paste(substr(as.character(ref),s-10,s-1),substr(as.character(ref),s,e),substr(as.character(ref),e+1,e+10),sep=" ")
    obs.tab <- count.cods(aln,s,e,res)
    cod.tab[j,] <- as.numeric(obs.tab[names(PHIX_GENETIC_CODE)])
    print(paste("done with ",tab$id[j],sep=""))
  }
}
cod.tab$sample <- tab$id
row.names(cod.tab) <- cod.tab$sample
#write.csv(cod.tab, "results/cod_tab_042219.csv", quote = F)
#write.csv(cod.tab, "results/cod_tab_g11.csv", quote = F)
#write.csv(cod.tab, "results/cod_tab_g45g128.csv", quote = F)



####process reads from 080218
sam.tab <- read.xlsx("rawdata/SubmissionSheet_Amplicon_v3.0_080218_jtv.xls", sheetIndex = 1)
sam.tab <- sam.tab[88:104,1:14]
names(sam.tab) <- as.character(as.matrix(sam.tab[1,]))
sam.tab <- sam.tab[-1,]

fns <- sort(list.files("rawdata/illumina_reads/080218/bam/", full.names = TRUE))
fnFs <- fns[grepl("sorted.bam.bai", fns)]
sample.names <- sapply(strsplit(fnFs, ".aligned"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])

#g41 <- data.frame(sample.names[str_detect(sample.names, "g41_t")])
tab <- data.frame(sample.names)
names(tab) <- "id"
tab$id <- str_replace(tab$id,"-","_")
tab$id <- str_replace(tab$id,"_redo","-redo")
tab$time <- str_split(tab$id,"_",simplify = T)[,2]
tab$rep <- str_split(tab$id,"_",simplify = T)[,3]
tab$site <- str_split(tab$id,"_",simplify = T)[,1]
tab$site <- str_extract(tab$site,"\\d+")
tab$id <- str_replace(tab$id,"_R1","")
tab$id <- str_replace(tab$id,"-redo","_redo")
tab$rep <- 'rep1'


##could not find a good tool to measure codon frequencies over time.
#have to write myself
#####try with bam files
###code works with single mix (not multiple sites mixed together)
cod.tab <- data.frame(matrix(nrow=nrow(tab),ncol=64))
names(cod.tab) <- names(PHIX_GENETIC_CODE)
slow <- F
if(slow){
  for(j in 1:nrow(tab)){
    tab$id[j]
    bam.sorted <- BamFile(paste("rawdata/illumina_reads/080218/bam/",tab$id[j],".aligned.sorted.bam",sep=""))
    #seqinfo(bam.sorted)
    #names(bam.sorted)
    #quickBamFlagSummary(bam.sorted)
    aln <- scanBam(bam.sorted)
    aln <- aln[[1]]
    #names(aln)
    #lapply(aln, function(xx) xx[1])
    cigs <- table(aln$cigar)
    cigs[cigs > length(aln$qname)*0.001]
    #filter alinged reads by
    #1 mapq
    #2 isize
    res <- as.numeric(tab$site[j])
    s <- res*3-2
    e <- res*3
    paste(substr(as.character(ref),s-10,s-1),substr(as.character(ref),s,e),substr(as.character(ref),e+1,e+10),sep=" ")
    obs.tab <- count.cods(aln,s,e,res)
    cod.tab[j,] <- as.numeric(obs.tab[names(PHIX_GENETIC_CODE)])
    print(paste("done with ",tab$id[j],sep=""))
  }
}
cod.tab$sample <- tab$id
row.names(cod.tab) <- cod.tab$sample
#write.csv(cod.tab, "results/cod_tab_080218.csv", quote = F)


cod.tab.1 <- read.csv("results/cod_tab_042219.csv")
cod.tab.2 <- read.csv("results/cod_tab_080218.csv")

cod.tab <- rbind(cod.tab.1, cod.tab.2)
#write.csv(cod.tab, "results/cod_tab_combined.csv", quote = F, row.names = F)



##naming was really bad. "pre" is at the start of the growth assay. 
##t0 is just before transformation. I double checked this in my lab notebook.


cod.tab.pre <- cod.tab[str_detect(cod.tab$X, "pre"),]
cod.tab.pre[is.na(cod.tab.pre)] <- 0 


##try dms_tools
##had to fix reads first
##pasted barcodes onto reads
#/mnt/c/Users/jt_laptop/Documents/perl/onetimers/add_barcodes.pl ../phix_41_43_106_R2.fastq rev > R2_fix.fastq
#/mnt/c/Users/jt_laptop/Documents/perl/onetimers/add_barcodes.pl ../phix_41_43_106_R1.fastq fwd > R1_fix.fastq
#mv R2_fix.fastq phix_fix_R2.fastq
##test on one sample before running batch

#dms2_bcsubamp --refseq ../../../refseq/G_wt.fa --alignspecs 1,528,12,35 --bclen 8 --maxmuts 7 --minq 20 --minreads 10 --minfraccall 0.90 --name bcsubamp --R1 phix_fix_R1.fastq --R2 phix_fix_R2.fastq

#dms2_diffsel --outdir dms2_diffsel --ncpus 10 --indir split --excludestop no --name g106_t0_rep1_dms2 

#dms2_batch_diffsel --outdir dms2_diffsel --ncpus 10 --indir split --excludestop no --batchfile 
