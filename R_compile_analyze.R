##code to pull data from various sources and join it

library(ggplot2)
library(stringr)
library(cowplot)
library(Biostrings)
library(seqinr)
library(reshape2)
library(dplyr)
library(performace)

PHIX_GENETIC_CODE <- GENETIC_CODE
attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"
ref <- readDNAStringSet("rawdata/G_wt.fasta")
refaa <- Biostrings::translate(ref)

#functions
{
  #function joins plaque size, sequencing, and viability data into one table
  merge_data <- function(site){
  plqs <- read.csv(file = paste("rawdata/plate_images/", site,"/", site, "_my_results_summary.csv",sep=""), row.names = 1, header = T)
  site2 <- tolower(site)
  seqs <- read.csv(file=paste("rawdata/sanger_reads/", site2, ".csv",sep=""), row.names=1,header=T)
  ####change shit here by hand if needed to make Lu happy
  ###eg
  ##seqs[seqs$sample.names=="G15-142",]$codon <- "TCC"
  ##seqs[seqs$sample.names=="G15-142",]$aa <- as.character(Biostrings::translate(DNAString(seqs[seqs$sample.names=="G15-142",]$codon)))
  foldx <- read.csv("results/ddG_foldX_GProtein.csv", header=T)
  foldxdimer <- read.csv("results/ddg_foldX_GProtein_dimer.csv", header=F)
  foldxall <- read.csv("results/ddG_foldX_GProtein.csv", sep="\t", header=T)
  site <- str_replace(site,"new","")
  res <- as.numeric(as.character(str_replace(site,"G","")))
  WTaa <- unname(as.character(subseq(refaa,res,res)))
  #reformat a bit
  plqs$sample.names <- str_split(plqs$plq,"_",simplify = T)[,1]
  plqs$sample.names <- str_replace(plqs$sample.names, "G", "")
  plqs$sample.names <- str_replace(plqs$sample.names, "new", "")
  seqs$sample.names <- str_replace(seqs$sample.names, "G", "")
  seqs$sample.names <- str_replace(seqs$sample.names, "new", "")
  
  #join sequencing frequency data, plaque size data, and foldX data into a table
  ##problem here with redo (new primers) samples because we did not redo plaq pics for previously observed codons
  join <- merge(seqs,plqs,by="sample.names") 
  seqs.short <- na.omit(seqs[!seqs$codon=="het",])
  seqs.short$codon <- factor(seqs.short$codon)
  aafreqs <- table(seqs.short$aa)/nrow(seqs.short)
  codfreq <- table(seqs.short$codon)/nrow(seqs.short)
  tablefreqs <- table(PHIX_GENETIC_CODE)/sum(table(PHIX_GENETIC_CODE))
  
  ##problem here with redo (new primers) samples because we did not redo plaq pics for previously observed codons
  #tmptab <- table(seqs.short[!seqs.short$codon %in% join$codon,]$codon)
  #add_cods <- names(tmptab[which(tmptab > 0)])
  #add <- data.frame(sample.names=4, plate.pos=NA, codon=add_cods, aa=PHIX_GENETIC_CODE[add_cods], plq=NA, radius.mean=NA, radius.median=NA, radius.getmode=NA, radius.getquant=NA, radius.max=NA, radius.min=NA, radius.st.err=NA, radius.length=NA)
  #join <- rbind(join, add)
  
  
  #join$seq
  join$aaexp <- NA
  join$aacount <-NA
  join$aafreq <- NA
  join$codonfreq <- NA
  join$fold <- NA
  join$bind <- NA
  
  for(i in 1:nrow(join)){
    print(i)
    join[i,]$aacount <- aafreqs[join[i,]$aa] * nrow(seqs.short)
    join[i,]$aafreq <- aafreqs[join[i,]$aa]
    join[i,]$codonfreq <- codfreq[join[i,]$codon]
    join[i,]$aaexp <- tablefreqs[join[i,]$aa] * nrow(seqs.short)
    get <- paste(WTaa,str_split(join[i,]$sample.names,"-",simplify=T)[,1],join[i,]$aa,sep="")
    if(join[i,]$aa == WTaa){
      join[i,]$fold <- 0
      join[i,]$bind <- 0
    }else{
      join[i,]$fold <- foldx[str_detect(foldx$mut,get),]$fold
      join[i,]$bind <- foldx[str_detect(foldx$mut,get),]$pentamer   #######Change to pentamer, trimer, dimer
    }
  }
  join$aanorm <- (sum(join$aacount)/sum(join$aacount/join$aaexp))*(join$aacount/join$aaexp)
  join_short <- join[,names(join) %in% c("sample.names", "codon", "aa", "radius.getquant", "aafreq", "aacount", "aaexp", "aanorm","codonfreq", "fold", "bind")]
  
  ##get fold/bind data from non-recovered variants
  foldx_short <- foldx[str_sub(foldx$mut,end=nchar(as.character(foldx$mut))-1) == paste(WTaa,res,sep=""),]
  foldx_short$aa <- str_sub(foldx_short$mut,-1,-1)
  foldx_short$sample.names <- res
  foldx_short$codon <- "XXX"
  #foldx_short$radius.getquant <- min(join_short$radius.getquant)
  foldx_short$radius.getquant <- NA
  foldx_short$aafreq <- NA 
  foldx_short$aacount <- NA 
  foldx_short$aaexp <- NA
  foldx_short$aanorm <- NA
  foldx_short$codonfreq <- 0
  #foldx_short$fold <- foldx_short$fold
  foldx_short$bind <- foldx_short$pentamer   #######Change to pentamer, trimer, dimer
  foldx_short <- foldx_short[,names(foldx_short) %in% names(join_short)]
  foldx_short <- foldx_short[,names(join_short)]
  foldx_short <- foldx_short[! foldx_short$aa %in% as.character(join_short$aa),]
  
  plot <- rbind(join_short,foldx_short)
  return(plot)
}
}

##compile plaque size measurements
{
dirs <- list.dirs("rawdata/plate_images/")
dirs_singles <- dirs[str_detect(dirs, "singles")]
dirs <- dirs[!str_detect(dirs, "singles")]
dirs <- dirs[!str_detect(dirs, "NNK")]
plaques <- unique(str_split(dirs,"\\/", simplify = T)[,4])
plaques <- plaques[-1]
plq.dat <- data.frame(matrix(ncol=9))
for(i in 1:length(plaques)){
  dat <- read.csv(paste("rawdata/plate_images/",plaques[i], "/", plaques[i], "_my_results_summary.csv",sep=""), row.names = 1)
  names(plq.dat) <- names(dat)
  plq.dat <- rbind(plq.dat, dat)
}
plq.dat <- na.omit(plq.dat)
}
#write.csv(plq.dat, "results/plaque_sizes.csv", quote = F, row.names = F)

#put together sequence counts, plaque sizes, etc. into one table (plot.mat.csv)
reanalyze <- FALSE
if(reanalyze) {
  g3plot <- merge_data("G3")
  g4plot <- merge_data("G4")
  g8plot <- merge_data("G8")
  g11plot <- merge_data("G11")
  g15plot <- merge_data("G15")
  #for 15-142, I called it a het, but Lu called it TCC.
  #had to go fix by hand.
  ##seqs[seqs$sample.names=="G15-142",]$codon <- "TCC"
  ##seqs[seqs$sample.names=="G15-142",]$aa <- as.character(Biostrings::translate(DNAString(seqs[seqs$sample.names=="G15-142",]$codon)))
  #g15plot <- plot
  g41plot <- merge_data("G41")
  #for 41-136, I failed the sequence because it is noisy, but Lu called it GCG.
  #had to go fix by hand.
  #site <- "G41"
  ##seqs[seqs$sample.names=="G41-136",]$codon <- "GCG"
  ##seqs[seqs$sample.names=="G41-136",]$aa <- as.character(Biostrings::translate(DNAString(seqs[seqs$sample.names=="G41-136",]$codon)))
  ##g41plot <- plot
  g43plot <- merge_data("G43")
  g45plot <- merge_data("G45")
  g72plot <- merge_data("G72")
  #oh man...have to check out all of these
  #c("72-11","72-22","72-3","72-58")
  #3 and 58 are probably fine, the other two have quite sketchy traces
  #drop them for now
  #site <- "G72"
  #join <- na.omit(join)
  #g72plot <- plot
  g74plot <- merge_data("G74")
  g80plot <- merge_data("G80")
  g81plot <- merge_data("G81")
  g88plot <- merge_data("G88")
  g106plot <-merge_data("G106")
  #site <- "G106"
  #106-26, 106-94 are probably just fine, but a little something strange with sequences.
  #drop them for now...i"m in a hurry
  #join <- na.omit(join)
  #g106plot <- plot
  g117plot <- merge_data("G117")
  g119plot <- merge_data("G119")
  #site <- "G119"
  #join <- na.omit(join)
  #g119plot <- plot
  g123plot <- merge_data("G123")
  g125plot <- merge_data("G125")
  g128plot <- merge_data("G128")
  g129plot <- merge_data("G129")
  g148plot <- merge_data("G148")
  ##g4plot_new <- must go make by hand
  ##g125plot_new <- must go make by hand
  ggplot(g11plot, aes(bind,fold, color=aafreq), label=aa) + 
    geom_point(aes(size=radius.getquant)) + 
    theme_minimal() +
    scale_color_gradient(low="#0091ff", high = "#f0650e", na.value = "black") +
    geom_text(aes(label=aa), hjust=1.5, vjust=1.5)

  plot.mat <- rbind(g3plot,
                    g4plot,
                    g8plot,
                    g11plot,
                    g15plot,
                    g41plot,
                    g43plot,
                    g45plot,
                    g72plot,
                    g74plot,
                    g80plot,
                    g81plot,
                    g88plot,
                    g106plot,
                    g117plot,
                    g119plot,
                    g123plot,
                    g125plot,
                    g128plot,
                    g129plot,
                    g148plot)
  
  plot.mat$site <- str_split(plot.mat$sample.names,"-",simplify = T)[,1]
  plot.mat$site <- as.factor(plot.mat$site)
  plot.mat$fill <- plot.mat$aafreq
  plot.mat$fill[!is.na(plot.mat$fill)] <- 1
  plot.mat[which(!is.na(plot.mat$aafreq)),]$fill <- as.character(plot.mat[which(!is.na(plot.mat$aafreq)),]$site)
  plot.mat$fill <- as.factor(plot.mat$fill)
  plot.mat[which(is.na(plot.mat$aafreq)),]$radius.getquant <- 0.1
  plot.mat$via <- as.factor(as.character(is.na(plot.mat$aafreq)))
  #write.csv(plot.mat,file="results/plot.mat.csv", quote=F, row.names = F)
  #write.csv(plot.mat,file="results/plot.mat_g4new.csv", quote=F, row.names = F)
  #write.csv(plot.mat,file="results/plot.mat_g125new.csv", quote=F, row.names = F)
}

#add to data a bit
plot.mat <- read.csv("results/plot.mat.csv")
plot.mat.g4new <- read.csv("results/plot.mat_g4new.csv", sep='\t')
plot.mat.g125new <- read.csv("results/plot.mat_g125new.csv")
start.freqs <- read.csv("results/start_freq_codon_combined.csv", row.names=1)
start.freqs.l <- read.csv("results/start_freq_codon_long.csv", row.names=1)
start.freqs.aa <- read.csv("results/start_freq_aa.csv", row.names=1)
row.names(start.freqs.l) <- paste(start.freqs.l$site,start.freqs.l$variable,sep="_")
row.names(start.freqs.aa) <- start.freqs.aa$Group.1
##remove old g4 and g125 and replace. must get plaque sizes from old ones. 
plot.new <- plot.mat[!plot.mat$site %in% c(4,125),]

plot.old <- plot.mat[plot.mat$site == 4,]
getcods <- plot.mat.g4new[is.na(plot.mat.g4new$radius.getquant),]$codon
#getcods <- getcods[!getcods %in% plot.old$codon]
lookup <- plot.old[plot.old$codon %in% getcods,]$radius.getquant
names(lookup) <- plot.old[plot.old$codon %in% getcods,]$codon
plot.mat.g4new[is.na(plot.mat.g4new$radius.getquant),]$radius.getquant <- lookup[plot.mat.g4new[is.na(plot.mat.g4new$radius.getquant),]$codon]

plot.old <- plot.mat[plot.mat$site == 125,]
getcods <- plot.mat.g125new[str_detect(plot.mat.g125new$sample.names,"-"),]$codon
#getcods <- getcods[!getcods %in% plot.old$codon]
lookup <- plot.old[plot.old$codon %in% getcods,]$radius.getquant
names(lookup) <- plot.old[plot.old$codon %in% getcods,]$codon
plot.mat.g125new[is.na(plot.mat.g125new$radius.getquant),]$radius.getquant <- lookup[plot.mat.g125new[is.na(plot.mat.g125new$radius.getquant),]$codon]

plot.mat <- rbind(plot.new, plot.mat.g125new, plot.mat.g4new)

plot.mat$match <- paste(plot.mat$site, plot.mat$codon, sep="_")
plot.mat$match2 <- paste(plot.mat$site, plot.mat$aa, sep="_")
plot.mat$start.cod <- start.freqs.l[plot.mat$match,]$value
plot.mat$start.aa <- start.freqs.aa[plot.mat$match2,]$x
cod.tab.counts <- table(PHIX_GENETIC_CODE)
plot.mat$start.aa.norm <- NA
plot.mat$aafreq.norm <- NA
for(i in 1:nrow(plot.mat)){
  plot.mat[i,]$start.aa.norm <- plot.mat[i,]$start.aa / cod.tab.counts[plot.mat[i,]$aa]
  plot.mat[i,]$aafreq.norm <- plot.mat[i,]$aafreq / cod.tab.counts[plot.mat[i,]$aa]
}
plot.mat[is.na(plot.mat$aafreq),]$aafreq <- 0
plot.mat[is.na(plot.mat$aafreq.norm),]$aafreq.norm <- 0
##only partially fixed. should have all codons listed, not just amino acids

###load in Jordan's data
rosetta <- read.csv("results/Trimer_Docking_PyRosetta_Yang.csv", stringsAsFactors = F)
names(rosetta) <- c("mutation","trimer_ddg")
rosetta_ddg <- rosetta
rosetta_ddg$trimer_ddg <- rosetta_ddg$trimer_ddg + 503.11
rosetta_ddg <- rosetta_ddg[!str_detect(rosetta_ddg$mutation,"wildtype"),]
rosetta_ddg$aa <- a(str_to_title(tolower(str_split(rosetta_ddg$mutation,"_",simplify = T)[,1])))
rosetta_ddg$site <- str_split(rosetta_ddg$mutation,"_",simplify = T)[,2]
for (i in 153:399){
   rosetta_ddg[i,]$aa <- str_split(rosetta_ddg[i,]$mutation,"",simplify = T)[,nchar(rosetta_ddg[i,]$mutation)]
   rosetta_ddg[i,]$site <- str_extract(rosetta_ddg[i,]$mutation,"\\d+")
}
row.names(rosetta_ddg) <- paste(rosetta_ddg$site,rosetta_ddg$aa,sep="_")
 #plot.mat$match3 <- paste(toupper(aaa(plot.mat$aa)),plot.mat$site,sep="_")
rosetta_fold <- read.csv("results/Combined_pyR_FoldX_data_v2.csv", stringsAsFactors = F)
index <- data.frame(site=1:176, seq=str_split(as.character(refaa), "", simplify=T)[1,])
row.names(rosetta_fold) <- rosetta_fold$mutation
plot.mat$pyfold <- rosetta_fold[paste(index[plot.mat$site,]$seq, plot.mat$site, plot.mat$aa, sep=""),]$Py.fold
plot.mat$rosetta <- rosetta_ddg[plot.mat$match2,]$trimer_ddg
plot.mat[is.na(plot.mat$rosetta),]$rosetta <- 0
#write.csv(plot.mat, file="results/data_combined.csv", quote=F, row.names = F)


datCodons <- read.csv("results/data_combined.csv", stringsAsFactors = F)
start_freq_codon <- read.csv("results/start_freq_codon.csv", row.names=1, stringsAsFactors = F)
start_freq_codon <- start_freq_codon[! row.names(start_freq_codon) %in% c("g4_lig_mix", "g125_lig_mix"),]
###getting noise from low read count data
###add cuttoff
sites <- start_freq_codon$site
#start_freq_codon <- start_freq_codon %>% replace(.<10, 0)
start_freq_codon$site <- sites

##analyze difference between synonymous codons
##datCodons is missing codon starting frequency for unobserved codons. 
###The datCodons table comes from the merge_data function above, but this
###function does not output codon frequencies for unobserved codons....need to fix
datCodons.short <- datCodons[datCodons$site %in% start_freq_codon$site,]
tmp <- datCodons.short[!is.na(datCodons.short$aaexp),]
tmp2 <- datCodons.short[is.na(datCodons.short$aaexp),]

#add unobserved codons
for(i in 1:nrow(tmp)){
  aa <- tmp[i,]$aa
  site <- tmp[i,]$site
  get <- which(PHIX_GENETIC_CODE == aa)
  add <- names(get)
  add <- add[!add %in% tmp[tmp$site == site,]$codon]
  if(!isEmpty(add)){
    new <- tmp[rep(i,length(add)), ]
    new$codon <- add
    new$match <- paste(new$site,add,sep="_")
    new$sample.names <- as.character(site)
    new$radius.getquant <- NA
    new$codonfreq <- NA
    new$start.cod <- NA
    new$via <- TRUE
    tmp <- rbind(tmp, new)
  }
}

##add unobserved amino acids
for(i in 1:nrow(tmp2)){
  aa <- tmp2[i,]$aa
  get <- which(PHIX_GENETIC_CODE == aa)
  add <- names(get)
  new <- tmp2[rep(i,length(add)), ]
  new$codon <- add
  new$match <- paste(new$site,add,sep="_")
}

##add freqs
start_freq_norm <- start_freq_codon
start_freq_norm[,2:65] <- start_freq_codon[,2:65]/start_freq_codon$rowsum
row.names(start_freq_norm) <- start_freq_norm$site
for(i in 1:nrow(tmp)){
  site <- as.character(tmp[i,]$site)
  cod <- as.character(tmp[i,]$codon)
  tmp[i,"start.cod"] <- start_freq_norm[site,cod]
}

#something went wrong. some codons are represented more than once
table(table(tmp$match) > 1)
fix <- names(which(table(tmp$match) > 1))
dups <- tmp[tmp$match %in% fix,]
View(dups)
###all appear to be instances where multiple pictures were taken
##will average plaque size
dup.means <- dups %>%
  group_by(match) %>%
  summarise(mean = mean(radius.getquant))
drop <- which(tmp$sample.names %in% dups$sample.names[c(2,4,6,8,10,12)])
tmp <- tmp[-drop,]
#write.csv(tmp, "results/data_combined_allcods.csv", quote = F, row.names = F)



##################################################################################################################
##this is a god damn mess. data_combine_allcods.csv is missing some codons that were observed in the illumina sequencing.
##just start over

start_freq <- read.csv("results/start_freq_codon.csv", row.names = 1)
start_freq$id <- row.names(start_freq)
start_freq <- start_freq[,!names(start_freq) %in% c('mapped','rowsum','site')]
start_freq_long <- melt(start_freq, id.vars = c('id'))
start_freq_long$site <- str_extract(start_freq_long$id, '\\d+')
names(start_freq_long) <- c('id', 'codon', 'readcnt' ,'site')

start_freq_long <- start_freq_long %>%
  group_by(id) %>%
  #summarise(n=n()) %>%
  mutate(readfreq = readcnt/sum(readcnt))

#merge g4 and g125 redo and old. average frequency. sum reads.
start_freq_long_merge <- start_freq_long[,-1]
start_freq_long_merge <- start_freq_long_merge %>%
  group_by(codon,site) %>%
  mutate(readfreq= mean(readfreq)) %>%
  mutate(readcnt=sum(readcnt))
start_freq_long_merge <- start_freq_long_merge[!duplicated(start_freq_long_merge),]

all_cods <- read.csv("results/all_cods.csv")
all_cods$site <- str_extract(all_cods$isolate, "\\d+")

geno_counts <- all_cods %>%
  group_by(site, codon) %>%
  summarise(plqcnt=n()) %>%
  mutate(genofreq = plqcnt/sum(plqcnt))

geno_counts$id <- paste(geno_counts$site, geno_counts$codon, sep="_")
start_freq_long_merge$id <- paste(start_freq_long_merge$site, start_freq_long_merge$codon, sep="_")

comp_freqs <- merge(geno_counts, start_freq_long_merge, by='id', all=T)
comp_freqs$codon <- str_split(comp_freqs$id, "_", simplify = T)[,2]
comp_freqs$site <- str_split(comp_freqs$id, "_", simplify = T)[,1]
comp_freqs <- comp_freqs[,c('id', 'plqcnt', 'genofreq', 'readcnt', 'readfreq', 'codon', 'site')]
comp_freqs$aa <- PHIX_GENETIC_CODE[comp_freqs$codon]
comp_freqs[is.na(comp_freqs)] <- 0
comp_freqs$site <- factor(comp_freqs$site, levels = sort(as.numeric(unique(comp_freqs$site))))

ntobs <- which(comp_freqs$genofreq == 0) 
ntobs <- c(ntobs, which(comp_freqs$readfreq == 0))
ntobs <- as.numeric(names(table(ntobs)[table(ntobs) ==1]))

comp_freqs$codonOE <- comp_freqs$genofreq-comp_freqs$readfreq
comp_freqs[ntobs,]$codonOE <- NA   #####not observed in sequence data

#get wt codons
cods <- data.frame(site = 1:176, seq = substring(ref, seq(1,nchar(ref),3), seq(3,nchar(ref),3)))
comp_freqs$wt <- cods[comp_freqs$site,]$seq
comp_freqs$wt <- comp_freqs$wt == comp_freqs$codon
comp_freqs$codonOE_plot <- NA
comp_freqs[ntobs,]$codonOE_plot <- comp_freqs[ntobs,]$genofreq-comp_freqs[ntobs,]$readfreq
comp_freqs[which(comp_freqs$codonOE_plot > 0),]$codonOE_plot <- max(na.omit(comp_freqs$codonOE))
comp_freqs[which(comp_freqs$codonOE_plot < 0),]$codonOE_plot <- min(na.omit(comp_freqs$codonOE))
comp_freqs[which(comp_freqs$codonOE == 0),]$codonOE <- NA



rosetta <- read.csv("rawdata/molecular_modeling/Trimer_Docking_PyRosetta_Yang.csv", stringsAsFactors = F)
names(rosetta) <- c("mutation","trimer_ddg")
rosetta_ddg <- rosetta
rosetta_ddg$trimer_ddg <- rosetta_ddg$trimer_ddg + 503.11
rosetta_ddg <- rosetta_ddg[!str_detect(rosetta_ddg$mutation,"wildtype"),]
rosetta_ddg$aa <- a(str_to_title(tolower(str_split(rosetta_ddg$mutation,"_",simplify = T)[,1])))
rosetta_ddg$site <- str_split(rosetta_ddg$mutation,"_",simplify = T)[,2]
for (i in 153:399){   ####Jordan changed naming formats half way though. Must deal with this.
  rosetta_ddg[i,]$aa <- str_split(rosetta_ddg[i,]$mutation,"",simplify = T)[,nchar(rosetta_ddg[i,]$mutation)]
  rosetta_ddg[i,]$site <- str_extract(rosetta_ddg[i,]$mutation,"\\d+")
}
row.names(rosetta_ddg) <- paste(rosetta_ddg$site,rosetta_ddg$aa,sep="_")

rosettaDimer <- read.csv("rawdata/molecular_modeling/Dimer_PyRosetta_Yang.csv", stringsAsFactors = F)
rosettaDimer_ddg <- rosettaDimer
names(rosettaDimer_ddg) <- c("mutation","dimer_ddg")
rosettaDimer_ddg$dimer_ddg <- rosettaDimer_ddg$dimer_ddg + 331.68
rosettaDimer_ddg  <- rosettaDimer_ddg[!str_detect(rosettaDimer_ddg$mutation,"wildtype"),]
rosettaDimer_ddg$aa <- a(str_to_title(tolower(str_split(rosettaDimer_ddg$mutation,"_",simplify = T)[,1])))
rosettaDimer_ddg$site <- str_extract(rosettaDimer_ddg$mutation, '\\d+')
rosettaDimer_ddg[which(is.na(rosettaDimer_ddg$aa)),]$aa <- str_split(rosettaDimer_ddg[which(is.na(rosettaDimer_ddg$aa)),]$mutation, '\\d+', simplify = T)[,2]
row.names(rosettaDimer_ddg) <- paste(rosettaDimer_ddg$site, rosettaDimer_ddg$aa, sep="_")

rosetta_join <- merge(rosetta_ddg, rosettaDimer_ddg, by=0, all=TRUE)


