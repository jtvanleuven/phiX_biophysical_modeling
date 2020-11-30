library(stringr)
library(Biostrings)
library(xlsx)
library(Rsamtools)
library(Rmisc)

PHIX_GENETIC_CODE <- GENETIC_CODE
attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"


count.cods <- function(aln,s,e){
  PHIX_GENETIC_CODE <- GENETIC_CODE
  attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"
  obs <- NULL
  for(i in 1:length(aln$qname)){
    seq <- aln$seq[i]
    qual <- aln$qual[i]
    cig <- aln$cigar[i]
    if(!is.na(aln$mapq[i])){    #some reads have no mapping quality
      if(aln$mapq[i] > 50){     #minimum mapping quality
        cig2 <- str_split(cig,"[A-Z]",simplify = T)
        cig3 <- str_split(cig,"\\d+",simplify = T)[,-1]
        if(s > aln$pos[i] && s < aln$pos[i]+sum(as.numeric(cig2[which(cig3=="M")]))){  ##make sure the codon of interest is within the mapping range
          if(all(!str_detect(cig, c("I","D","P")))){   #ignore reads with indels
            if(cig3[1] == "S" | cig3[1] == "H"){
              s2 <- (s - aln$pos[i]) + as.numeric(cig2[1]) + 1
            #}else if(cig3[1] == "I"){
              #s2 <- (s - aln$pos[i]) + as.numeric(cig2[1]) + 1
            }else if(cig3[1] == "M"){
              s2 <- s - (aln$pos[i] - 1)
            }else{
              print(paste("unreconginzed cigar:", cig, sep=" "))
            }
            e2 <- s2+2
            cod <- substr(as.character(seq),s2,e2)
            #paste(substr(as.character(seq),s2-10,s2-1),substr(as.character(seq),s2,e2),substr(as.character(seq),e2+1,e2+10),sep=" ")
            cod.qual <- substr(as.character(qual),s2,e2)
            #paste(substr(as.character(qual),s2-10,s2-1),substr(as.character(qual),s2,e2),substr(as.character(qual),e2+1,e2+10),sep=" ")
            #as.numeric(as.character(charToRaw(cod.qual))) - 33 + 24 #not really sure why R doesn't get ascii right
            #qual[cod]
            #aln$pos[i]
            obs <- c(obs,cod)
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



##on linux for processing merged reads
#java -jar /mnt/c/Users/jt_laptop/Downloads/keep/Trimmomatic-0.38/trimmomatic-0.3jar PE -threads 10 phix_41_43_106_R1.fastq phix_41_43_106_R2.fastq phix_41_43_106_F_paired.fq phix_41_43_106_F_unpaired.fq phix_41_43_106_R_paired.fq phix_41_43_106_R_unpaired.fq ILLUMINACLIP:MiSeqIBEST-PE.fa:2:20:10 LEADING:3 TRAILING:3 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:50
#Input Read Pairs: 3238372 Both Surviving: 2865996 (88.50%) Forward Only Surviving: 318308 (9.83%) Reverse Only Surviving: 26931 (0.83%) Dropped: 27137 (0.84%)
#bwa index -a is G_wt.fa
#samtools faidx G_wt.fa
# bwa mem -t 10 ../refseq/G_wt.fa reads/phix_41_43_106_F_paired.fq reads/phix_41_43_106_R_paired.fq > aln-pe.sam
#samtools view -bS aln-pe.sam -o aln-pe.bam
#samtools sort aln-pe.bam aln-pe.sorted
#samtools mpileup -d 10000000 -f ../refseq/G_wt.fa aln-pe.sorted.bam > mpileup-pe.txt
#java -jar /mnt/c/Users/jt_laptop/Downloads/keep/VarScan.v2.3.9.jar  mpileup2snp-pe.txt --output-vcf mpileup-pe.varscan --strand-filter 0
#/mnt/c/Users/jt_laptop/Documents/perl/parse_pileup_v4.pl mpileup-pe.txt 250 > mpileup-pe.parsed


##on linux for demultiplexing reads
# module load python
# source /mnt/lfs2/grcuser/opt/dbcAmpVenv/bin/activate
# splitReadsBySample.py -1 phix_41_43_106_R1.fastq -2 phix_41_43_106_R2.fastq
# run_bwa.sh
# run_samtools.sh

ref <- readDNAStringSet("rawdata/G_wt.fasta")
sam.tab <- read.xlsx("rawdata/SubmissionSheet_Amplicon_v3.0_042219_jtv.xls", sheetIndex = 1)
sam.tab <- sam.tab[15:157,1:14]
names(sam.tab) <- as.character(as.matrix(sam.tab[1,]))
sam.tab <- sam.tab[-1,]

fns <- sort(list.files("rawdata/illumina_reads/042219/", full.names = TRUE))
fnFs <- fns[grepl("R1.fastq.gz", fns)]
pileup <- fns[grepl("parse", fns)]
sample.names <- sapply(strsplit(fnFs, ".fastq"), `[`, 1)
sample.names <- c(sapply(strsplit(sample.names, "/"), `[`, str_count(sample.names,"/")+1)[1,])

#g41 <- data.frame(sample.names[str_detect(sample.names, "g41_t")])
tab <- data.frame(sample.names)
names(tab) <- "id"
tab$id <- str_replace(tab$id,"redo_","")
tab$id <- str_replace(tab$id,"-","_")
tab$id <- str_replace(tab$id,"_R1","")
tab$time <- str_split(tab$id,"_",simplify = T)[,2]
tab$rep <- str_split(tab$id,"_",simplify = T)[,3]
tab$site <- str_split(tab$id,"_",simplify = T)[,1]
tab$site <- str_extract(tab$site,"\\d+")
tab$id <- str_replace(tab$id,"_R1","")
tab[which(str_detect(tab$id, "NTC")),]$time <- 90 
tab[which(str_detect(tab$id, "NTC")),]$rep <- c(2,3,1)


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
    bam.sorted <- BamFile(paste("rawdata/illumina_reads/042219/",tab$id[j],"_bwa_mapped.sorted.bam",sep=""))
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
    obs.tab <- count.cods(aln,s,e)
    cod.tab[j,] <- as.numeric(obs.tab[names(PHIX_GENETIC_CODE)])
    print(paste("done with ",tab$id[j],sep=""))
  }
}
cod.tab$sample <- tab$id
row.names(cod.tab) <- cod.tab$sample
write.csv(cod.tab, "cod_tab.csv", quote = F)
cod.tab <- read.csv("cod_tab.csv")

cod.tab.pre <- cod.tab[str_detect(cod.tab$X, "pre"),]
cod.tab.pre[is.na(cod.tab.pre)] <- 0 








setwd("G41_G43_G106_exp/")
#cod.tab <- readRDS("cod.tab.rds")
cod.tab <- cbind(tab,cod.tab)
#add mapping data
cod.tab$mapped <- NA
for(i in 1:nrow(cod.tab)){
  cod.tab$id[i]
  bam.sorted <- BamFile(paste("reads/split/",cod.tab$id[i],"_bwa.sorted.bam",sep=""))
  aln <- scanBam(bam.sorted)
  aln <- aln[[1]]
  cod.tab[i,]$mapped <- as.numeric(table(aln$mapq > 50)[2])
}
cod.tab[is.na(cod.tab)] <- 0
cod.tab$rowsum <- rowSums(cod.tab[,5:68])
#saveRDS(cod.tab,"cod.tab.new.split2.rds")
cod.tab.new <- readRDS("cod.tab.new.rds")
cod.tab.new.split2 <- readRDS("cod.tab.new.split2.rds")


cod.tab2 <- cod.tab[cod.tab$rowsum > 500,]              ##########################dropping samples here
cod.tab2 <- cod.tab
cod.tab2[,5:68] <- cod.tab2[,5:68]/cod.tab2$rowsum
#saveRDS(cod.tab2,"cod.tab2.rds")
cod.tab2 <- readRDS("cod.tab2.rds")


lig2 <- cod.tab2[str_detect(cod.tab2$id, "lig"),]
lig2 <- lig2[,c(-2,-3)]
lig <- cod.tab[str_detect(cod.tab$id, "lig"),]
lig <- lig[,c(-2,-3)]
#write.csv(lig,"start_freq_codon.csv")
lig2 <- lig2[,!colnames(lig2) %in% c("mapped","rowsum")]
lig.long <- melt(lig2,id=c("id","site"))
ggplot(lig.long, aes(x=site, y=value, fill=variable))+
  geom_bar(stat="identity")


ggplot(lig.long, aes(x=value)) +
  geom_vline(xintercept = 0.015625) +
  geom_histogram(bins = 20) +
  facet_wrap(~site, scales = "free") +
  theme_grey(base_size = 10)



###probably want to convert to amino acid

lig.long$aa <- PHIX_GENETIC_CODE[lig.long$variable]
lig.long$ag <- paste(lig.long$site,lig.long$aa,sep="_")
lig.long.aa <- aggregate(lig.long$value, by=list(lig.long$ag), FUN=sum)
lig.long.aa$site <- str_split(lig.long.aa$Group.1,"_",simplify = T)[,1]
lig.long.aa$aa <- str_split(lig.long.aa$Group.1,"_",simplify = T)[,2]
#write.csv(lig.long,file="start_freq_codon_long.csv")
#write.csv(lig.long.aa,file="start_freq_aa_long.csv")
#lig.long <- read.csv("start_freq_aa_long.csv")

ggplot(lig.long, aes(x=as.factor(site), y=value, group=variable, color=variable))+
  geom_line(position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=8), axis.title.x = element_blank()) +
  ylab("Codon frequency")
  #ylim(0,0.1)

ggplot(lig.long, aes(x=value, y=variable, group=site, color=site))+
  geom_line(position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=8), axis.title.x = element_blank()) +
  ylab("Codon frequency")
#ylim(0,0.1)


lig.long.aa <- lig.long.aa[order(lig.long.aa$site),]
ggplot(lig.long.aa, aes(x=x))+
  geom_vline(xintercept = 0.04761905) +
  geom_histogram(bins = 12) +
  facet_wrap(~site) +
  #facet_wrap(~site, scales = "free") +
  theme_grey(base_size = 10) +
  theme_bw() +
  xlab("Amino acid frequency")


ggplot(lig.long.aa, aes(x=site, y=x), label=aa) +
  geom_hline(yintercept =  0.04761905) +
  #ylim(0,0.11) +
  geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=2)



cod.tab3 <- cod.tab2[!str_detect(cod.tab2$id,"mix"),]

#site 41
{
g41 <- cod.tab3[str_detect(cod.tab3$id, "g41_t"),]
g41 <- g41[,c(-2,-3)]
g41 <- g41[,!colnames(g41) %in% c("mapped","rowsum")]
g41.long <- melt(g41,id=c("id","site"))
g41.long$aa <- PHIX_GENETIC_CODE[g41.long$variable]
g41.long$ag <- paste(g41.long$id,g41.long$aa,sep="_")
g41.long.aa <- aggregate(g41.long$value, by=list(g41.long$ag), FUN=sum)
g41.long.aa$site <- str_split(g41.long.aa$Group.1,"_",simplify = T)[,1]
g41.long.aa$aa <- str_split(g41.long.aa$Group.1,"_",simplify = T)[,4]
g41.long.aa$time <- str_split(g41.long.aa$Group.1,"_",simplify = T)[,2]
g41.long.aa$rep <- str_split(g41.long.aa$Group.1,"_",simplify = T)[,3]
g41.long.aa$time <- str_split(g41.long.aa$time,"",simplify = T)[,2]

g41.lig.long <- lig.long.aa[lig.long.aa$site=="41",]
g41.lig.long$time <- -1
g41.lig.long$rep="rep1"

g41.long.aa <- rbind(g41.long.aa,g41.lig.long)

ggplot(g41.long.aa, aes(x=time, y=x), label=aa) +
  geom_hline(yintercept =  0.04761905) +
  ylim(0,0.11) +
  geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=3)

g41_summary <- summarySE(g41.long.aa, measurevar="x", groupvars=c("aa","time"))
g41_summary[is.na(g41_summary)] <- 0

ggplot(g41_summary, aes(x=time, y=x, group=aa, color=aa)) +
         geom_line(position=position_dodge(0.2)) +
         geom_point(size=3, position=position_dodge(0.2)) +
         geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2))
}      


#site 43
{
  g43 <- cod.tab3[str_detect(cod.tab3$id, "g43_t"),]
  g43 <- g43[!str_detect(g43$id, "g43_t2_rep1"),]    ###I clearly mixed up the samples.....this is a G41
  g43 <- g43[,c(-2,-3)]
  g43 <- g43[,!colnames(g43) %in% c("mapped","rowsum")]
  g43.long <- melt(g43,id=c("id","site"))
  g43.long$aa <- PHIX_GENETIC_CODE[g43.long$variable]
  g43.long$ag <- paste(g43.long$id,g43.long$aa,sep="_")
  g43.long.aa <- aggregate(g43.long$value, by=list(g43.long$ag), FUN=sum)
  g43.long.aa$site <- str_split(g43.long.aa$Group.1,"_",simplify = T)[,1]
  g43.long.aa$aa <- str_split(g43.long.aa$Group.1,"_",simplify = T)[,4]
  g43.long.aa$time <- str_split(g43.long.aa$Group.1,"_",simplify = T)[,2]
  g43.long.aa$rep <- str_split(g43.long.aa$Group.1,"_",simplify = T)[,3]
  g43.long.aa$time <- str_split(g43.long.aa$time,"",simplify = T)[,2]
  
  g43.lig.long <- lig.long.aa[lig.long.aa$site=="43",]
  g43.lig.long$time <- -1
  g43.lig.long$rep="rep1"
  
  g43.long.aa <- rbind(g43.long.aa,g43.lig.long)
  
  ggplot(g43.long.aa, aes(x=time, y=x), label=aa) +
    geom_hline(yintercept =  0.015625) +
    geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4))
  
  ggplot(g43.long.aa, aes(x=time, y=x), label=aa) +
    geom_hline(yintercept =  0.04761905) +
    ylim(0,0.11) +
    geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=3)
  
  g43_summary <- summarySE(g43.long.aa, measurevar="x", groupvars=c("aa","time"))
  g43_summary[is.na(g43_summary)] <- 0
  
  numcols <- length(unique(g43_summary$aa))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  
  ggplot(g43_summary, aes(x=time, y=x, group=aa, color=aa)) +
    #geom_hline(yintercept =  0.04761905) +
    geom_line(position=position_dodge(0.2), size=1.2, alpha=0.8) +
    #geom_point(size=3, position=position_dodge(0.2), alpha=0.8) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
    #ylim(0,0.12) +
    ylab("Genotype frequency") +
    scale_fill_manual(values = getPalette(numcols)) +
    scale_color_manual(values = getPalette(numcols))

} 


#site 106
{
  g106 <- cod.tab3[str_detect(cod.tab3$id, "g106_t"),]
  g106 <- g106[!str_detect(g106$id, "g106_t2_rep1"),]    ###I clearly mixed up the samples.....this is a G41
  g106 <- g106[,c(-2,-3)]
  g106 <- g106[,!colnames(g106) %in% c("mapped","rowsum")]
  g106.long <- melt(g106,id=c("id","site"))
  g106.long$aa <- PHIX_GENETIC_CODE[g106.long$variable]
  g106.long$ag <- paste(g106.long$id,g106.long$aa,sep="_")
  g106.long.aa <- aggregate(g106.long$value, by=list(g106.long$ag), FUN=sum)
  g106.long.aa$site <- str_split(g106.long.aa$Group.1,"_",simplify = T)[,1]
  g106.long.aa$aa <- str_split(g106.long.aa$Group.1,"_",simplify = T)[,4]
  g106.long.aa$time <- str_split(g106.long.aa$Group.1,"_",simplify = T)[,2]
  g106.long.aa$rep <- str_split(g106.long.aa$Group.1,"_",simplify = T)[,3]
  g106.long.aa$time <- str_split(g106.long.aa$time,"",simplify = T)[,2]
  
  g106.lig.long <- lig.long.aa[lig.long.aa$site=="106",]
  g106.lig.long$time <- -1
  g106.lig.long$rep="rep1"
  
  g106.long.aa <- rbind(g106.long.aa,g106.lig.long)
  
  ggplot(g106.long.aa, aes(x=time, y=x), label=aa) +
    geom_hline(yintercept =  0.015625) +
    geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4))
  
  ggplot(g106.long.aa, aes(x=time, y=x), label=aa) +
    geom_hline(yintercept =  0.04761905) +
    ylim(0,0.11) +
    geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=3)
  
  g106_summary <- summarySE(g106.long.aa, measurevar="x", groupvars=c("aa","time"))
  g106_summary[is.na(g106_summary)] <- 0
  
  ggplot(g106_summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_hline(yintercept =  0.04761905) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
    #ylim(0,0.12) +
    ylab("Genotype frequency")
} 


##put together sites to get a really rough idea of fitness from freq change
all <- rbind(g106.long.aa, g41.long.aa, g43.long.aa)
all <- all[all$time == c(1,5),]
all$Group.1 <- str_replace(all$Group.1,"t1_","")
all$Group.1 <- str_replace(all$Group.1,"t5_","")
all.w <- dcast(all, Group.1 + site + aa + rep ~ time, value.var = "x")
#something messed up with site 41
all.w <- all.w[!all.w$site=="g41",]
all.w$fit <- (all.w$`5`-all.w$`1`)/all.w$`1`
hist(all.w$fit,breaks=20)
write.csv(all.w,"fit_freq_g43_106.csv")

###should work the same, so long as it knows what site to go get
###just loop through the same data 3 times...one for each site
{
mix.41 <- mix
mix.41$site <- "41"
mix.43 <- mix
mix.43$site <- "43"
mix.106 <- mix
mix.106$site <- "106"


tab <- rbind(mix.41,mix.43,mix.106)

cod.tab <- data.frame(matrix(nrow=nrow(tab),ncol=64))
names(cod.tab) <- names(PHIX_GENETIC_CODE)
for(j in 1:nrow(tab)){
  tab$id[j]
  bam.sorted <- BamFile(paste("reads/split2/",tab$id[j],"_bwa.sorted.bam",sep=""))
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
  obs.tab <- count.cods(aln,s,e)
  cod.tab[j,] <- as.numeric(obs.tab[names(PHIX_GENETIC_CODE)])
  print(paste("done with ",tab$id[j],sep=""))
}

cod.tab$mapped <- NA
cod.tab <- cbind(tab,cod.tab)
for(i in 1:nrow(cod.tab)){
  cod.tab$id[i]
  bam.sorted <- BamFile(paste("reads/split/",codtab$id[i],"_bwa.sorted.bam",sep=""))
  aln <- scanBam(bam.sorted)
  aln <- aln[[1]]
  cod.tab[i,]$mapped <- as.numeric(table(aln$mapq > 50)[2])
}

cod.tab[is.na(cod.tab)] <- 0
cod.tab$rowsum <- rowSums(cod.tab[,5:68])
#saveRDS(cod.tab,"cod.tab.mix.rds")
cod.tab.mix <- readRDS("cod.tab.mix.rds")

cod.tab2.mix <- cod.tab.mix[cod.tab.mix$rowsum > 500,]              ##########################dropping samples here
#cod.tab2 <- cod.tab
cod.tab2.mix[,5:68] <- cod.tab2.mix[,5:68]/cod.tab2.mix$rowsum
#saveRDS(cod.tab2.mix,"cod.tab2.mix.rds")
cod.tab2.mix <- readRDS("cod.tab2.mix.rds")


mix2 <- cod.tab2.mix[str_detect(cod.tab2.mix$id, "mix"),]
mix2 <- mix2[,c(-2,-3)]
mix <- cod.tab.mix[str_detect(cod.tab.mix$id, "mix"),]
mix <- mix[,c(-2,-3)]
#write.csv(lig,"start_freq_codon.csv")
mix2 <- mix2[,!colnames(mix2) %in% c("mapped","rowsum")]
mix.long <- melt(mix2,id=c("id","site"))


#compare mix to singles.....somehow
mypalette <- c('#4d4d4d','#878787','#8dd3c7','#003c30','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
pie(rep(1,length(mypalette)), col=mypalette, labels=mypalette)
##start with g41
{
cod.tab3.mix <- cod.tab2.mix[cod.tab2.mix$site == "41",]
g41.mix <- cod.tab3.mix[,c(-2,-3)]
g41.mix <- g41.mix[,!colnames(g41.mix) %in% c("mapped","rowsum")]
g41.mix.long <- melt(g41.mix,id=c("id","site"))
g41.mix.long$aa <- PHIX_GENETIC_CODE[g41.mix.long$variable]
g41.mix.long$ag <- paste(g41.mix.long$id,g41.mix.long$aa,sep="_")
g41.mix.long.aa <- aggregate(g41.mix.long$value, by=list(g41.mix.long$ag), FUN=sum)
g41.mix.long.aa$site <- "41"
g41.mix.long.aa$aa <- str_split(g41.mix.long.aa$Group.1,"_",simplify = T)[,4]
g41.mix.long.aa$time <- str_split(g41.mix.long.aa$Group.1,"_",simplify = T)[,2]
g41.mix.long.aa$rep <- str_split(g41.mix.long.aa$Group.1,"_",simplify = T)[,3]
g41.mix.long.aa$time <- str_split(g41.mix.long.aa$time,"",simplify = T)[,2]

g41.mix.long.aa <- rbind(g41.mix.long.aa,g41.lig.long)

ggplot(g41.mix.long.aa, aes(x=time, y=x), label=aa) +
  geom_hline(yintercept =  0.04761905) +
  #ylim(0,0.11) +
  geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=3)

g41.mix.summary <- summarySE(g41.mix.long.aa, measurevar="x", groupvars=c("aa","time"))
g41.mix.summary[is.na(g41.mix.summary)] <- 0

ggplot(g41.mix.summary, aes(x=time, y=x, group=aa, color=aa)) +
  geom_line(position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2))


p1_g41 <- ggplot(g41_summary, aes(x=time, y=x, group=aa, color=aa)) +
  geom_line(position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
  ylab("") +
  xlab("")

  
p2_g41 <- ggplot(g41.mix.summary, aes(x=time, y=x, group=aa, color=aa)) +
  geom_line(position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
  ylab("") +
  xlab("")

legend <- get_legend(p1_g41 + theme(legend.position="bottom"))

prow <- plot_grid(
  p1_g41+theme(legend.position="none"),
  p2_g41+theme(legend.position="none"),
  labels = c("A", "B"))

plot_grid(prow,legend,ncol=1,rel_heights = c(1,0.2))
}

##g43
{
  cod.tab3.mix <- cod.tab2.mix[cod.tab2.mix$site == "43",]
  g43.mix <- cod.tab3.mix[,c(-2,-3)]
  g43.mix <- g43.mix[,!colnames(g43.mix) %in% c("mapped","rowsum")]
  g43.mix.long <- melt(g43.mix,id=c("id","site"))
  g43.mix.long$aa <- PHIX_GENETIC_CODE[g43.mix.long$variable]
  g43.mix.long$ag <- paste(g43.mix.long$id,g43.mix.long$aa,sep="_")
  g43.mix.long.aa <- aggregate(g43.mix.long$value, by=list(g43.mix.long$ag), FUN=sum)
  g43.mix.long.aa$site <- "43"
  g43.mix.long.aa$aa <- str_split(g43.mix.long.aa$Group.1,"_",simplify = T)[,4]
  g43.mix.long.aa$time <- str_split(g43.mix.long.aa$Group.1,"_",simplify = T)[,2]
  g43.mix.long.aa$rep <- str_split(g43.mix.long.aa$Group.1,"_",simplify = T)[,3]
  g43.mix.long.aa$time <- str_split(g43.mix.long.aa$time,"",simplify = T)[,2]
  
  g43.mix.long.aa <- rbind(g43.mix.long.aa,g43.lig.long)
  
  ggplot(g43.mix.long.aa, aes(x=time, y=x), label=aa) +
    geom_hline(yintercept =  0.04761905) +
    #ylim(0,0.11) +
    geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=3)
  
  g43.mix.summary <- summarySE(g43.mix.long.aa, measurevar="x", groupvars=c("aa","time"))
  g43.mix.summary[is.na(g43.mix.summary)] <- 0
  
  ggplot(g43.mix.summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2))
  
  
  p1_g43 <- ggplot(g43_summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
    ylab("") +
    xlab("")
  
  
  p2_g43 <- ggplot(g43.mix.summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
    ylab("") +
    xlab("")
  
  legend <- get_legend(p1_g43 + theme(legend.position="bottom"))
  
  prow <- plot_grid(
    p1_g43+theme(legend.position="none"),
    p2_g43+theme(legend.position="none"),
    labels = c("A", "B"))
  
  plot_grid(prow,legend,ncol=1,rel_heights = c(1,0.2))
}


##g106
{
  cod.tab3.mix <- cod.tab2.mix[cod.tab2.mix$site == "106",]
  g106.mix <- cod.tab3.mix[,c(-2,-3)]
  g106.mix <- g106.mix[,!colnames(g106.mix) %in% c("mapped","rowsum")]
  g106.mix.long <- melt(g106.mix,id=c("id","site"))
  g106.mix.long$aa <- PHIX_GENETIC_CODE[g106.mix.long$variable]
  g106.mix.long$ag <- paste(g106.mix.long$id,g106.mix.long$aa,sep="_")
  g106.mix.long.aa <- aggregate(g106.mix.long$value, by=list(g106.mix.long$ag), FUN=sum)
  g106.mix.long.aa$site <- "106"
  g106.mix.long.aa$aa <- str_split(g106.mix.long.aa$Group.1,"_",simplify = T)[,4]
  g106.mix.long.aa$time <- str_split(g106.mix.long.aa$Group.1,"_",simplify = T)[,2]
  g106.mix.long.aa$rep <- str_split(g106.mix.long.aa$Group.1,"_",simplify = T)[,3]
  g106.mix.long.aa$time <- str_split(g106.mix.long.aa$time,"",simplify = T)[,2]
  
  g106.mix.long.aa <- rbind(g106.mix.long.aa,g106.lig.long)
  
  ggplot(g106.mix.long.aa, aes(x=time, y=x), label=aa) +
    geom_hline(yintercept =  0.04761905) +
    #ylim(0,0.11) +
    geom_jitter(aes(color=aa), position = position_jitterdodge(dodge.width=0.1, jitter.width=0.4), size=3)
  
  g106.mix.summary <- summarySE(g106.mix.long.aa, measurevar="x", groupvars=c("aa","time"))
  g106.mix.summary[is.na(g106.mix.summary)] <- 0
  
  ggplot(g106.mix.summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2))
  
  
  p1_g106 <- ggplot(g106_summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
    ylab("") +
    xlab("")
  
  
  p2_g106 <- ggplot(g106.mix.summary, aes(x=time, y=x, group=aa, color=aa)) +
    geom_line(position=position_dodge(0.2)) +
    geom_point(size=3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin=x-sd, ymax=x+sd), position=position_dodge(0.2)) +
    ylab("") +
    xlab("")
    #ylim(c(0,0.1))
  
  legend <- get_legend(p1_g106 + theme(legend.position="bottom"))
  
  plot_grid(
    p1_g106+theme(legend.position="none"),
    p2_g106+theme(legend.position="none"),
    p2_g106+theme(legend.position="none") + ylim(c(0,0.06)),
    legend,
    labels = c("A", "B", "C"))
  
  #plot_grid(prow,legend,nc,rel_heights = c(1,0.2))
}

}







#parse mpileup
# for (i in 1:nrow(sing)){
#   file <- read.csv(pileup[str_detect(pileup, sing[i,]$id)], header = 0)
#   names(file) <- c("numvars","count2","sum","nA","nT","nG","nC","na","nt","ng","nc","nmatch","nopp")
#   
# }
# names(sing) <- c("id","time","rep","site","covA","covB","covC")


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
