library(ggplot2)
library(stringr)
library(cowplot)
library(Biostrings)
library(seqinr)
library(reshape2)
library(dplyr)
library(viridis)


#functions
{
  in.seq <- function(x) {
    # returns TRUE for elments within ascending sequences
    (c(diff(x, 1), NA) == 1 & c(NA, diff(x,2), NA) == 2)
  tmp}
  
  contractSeqs <-  function(x) {
    # returns string formatted with contracted sequences
    x[in.seq(x)] <- ""
    gsub(",{2,}", "-", paste(x, collapse=","), perl=TRUE)
  }
  
  findIntRuns <- function(run){
    rundiff <- c(1, diff(run))
    difflist <- split(run, cumsum(rundiff!=1))
    unlist(lapply(difflist, function(x){
      if(length(x) %in% 0:1) as.character(x) else paste0(x[1], "-", x[length(x)])
    }), use.names=FALSE)
  }
}

ref <- readDNAStringSet("rawdata/G_wt.fasta")
refaa <- Biostrings::translate(ref)

##problems
###depends on plaque size code from old stuff still. currently reads in "plot.mat.csv" which is created in "R_visualize_data_v2.R"


plot.mat <- read.csv("results/plot.mat.csv", stringsAsFactors = F) 
foldx <- read.csv(file="results/ddG_foldX_GProtein.csv", header=T, stringsAsFactors = F)
foldx2 <- foldx
foldx2$aa <- str_sub(foldx2$mut,1,1)
foldx2$var <- str_sub(foldx2$mut,-1,-1)
foldx2$site <- str_extract(foldx2$mut,'\\d+')

##pretty sure craig made this table
dat <- read.csv("results/Combined_pyR_FoldX_data_v3.csv", row.names = 1, stringsAsFactors = F)
##pretty sure I made this table. Contains plaque size and codon recovery data.
datCodons <- read.csv("results/data_combined_allcods.csv", stringsAsFactors = F)

seq.counts <- read.csv(file="results/protein_G_mutagenesis_AA_counts.csv", header=T)
##F-G data
GFbind <- read.csv("results/G_mut_GF_bind_ddg_stdev.csv", header = F, stringsAsFactors = F)

dat <- read.csv("results/Combined_pyR_FoldX_data_v3.csv", row.names = 1, stringsAsFactors = F)


#H = α-helix
#B = residue in isolated β-bridge
#E = extended strand, participates in β ladder
#G = 3-helix (310 helix)
#I = 5 helix (π-helix)
#T = hydrogen bonded turn
#S = bend 
#solvent accessibility from yesol
##she ran dssp (Kabsch,W. and Sander,C. (1983) Biopolymers 22, 2577-2637)
dssp <- read.csv("results/phixG_dssp.csv", row.names = 1)
##don't really like the results from this. 
###residues with high values don't seem to have particularly large RSA values

##tried another surface calculator: http://cib.cf.ocha.ac.jp/bitool/ASA/
###used pentamer structure
rsa2 <- read.table("rawdata/surface_pentamer.tab", sep="\t", header = T)
names(rsa2) <- c("chain","site","aa","area","rsa")
rsa2 <- rsa2[rsa2$chain %in% c("F","G","H","I","J"),]  ##G chains
rsa_sum <- rsa2 %>%
  filter(chain %in% c("F","G","H","I","J")) %>%
  group_by(site) %>%
  summarize(mean=mean(rsa))
names(rsa_sum) <- c("RESIDUE", "RSA2")

tmp <- merge(dssp, rsa_sum, by="RESIDUE") 
ggplot(tmp, aes(x=RSA, y=RSA2)) +
  geom_hline(yintercept = .5) +
  geom_vline(xintercept = .5) +
  geom_point(size=2, alpha=0.5) +
  theme_classic()
##hmmm, maybe yesol used g-only protien



##sort of a pain to get rate4site to work
##needs installed on windows machine
##takes *.fasta alignment (not nexus)
##takes newick alignment saved from figtree
Grates <- read.table("results/G_rate2site.txt", header=F, row.names = 1, stringsAsFactors = F)
ggplot(Grates, aes(x=V3, y=V5)) + geom_point()
names(Grates) <- c("wtAA","rate4site", "err", "raw", "something")
##redid everything with more sequences. Used blastp to search for G proteins
##took all sequences with evalue < 1E-10
gatherseq <- read.fasta("rawdata/lcl_Query_51641_209_other_sequences.aln", seqtype = "AA", as.string = T)
##cut out sequences that are less than 75% in length
longs <- str_count(as.character(gatherseq),"[A-Z]") > 0.75*175
gatherseq.sub <- gatherseq[longs]
dups <- duplicated(as.character(gatherseq.sub))
gatherseq.sub <- gatherseq.sub[!dups]
#write.fasta(sequences=gatherseq.sub, names = names(gatherseq.sub), file.out = "results/blastsearch.fasta")

##simplified blastsearch.fasta sequences by collapsing 99% identical sequencing and removing sequences that are not full length
#realigned using Muscle v3.8
#outputed conservation score from jalview
consv <- read.csv("results/jalview_consv_99ident_align", header=F, row.names = 1, stringsAsFactors = F)
consv <- consv[-176] ##stop is NA
Grates$jalview <- t(consv)
Grates$site <- 1:175
ggplot(Grates, aes(x=rate4site, y=jalview)) + geom_point()
ggplot(Grates, aes(x=raw, y=jalview)) + geom_point()
plot.cons <- Grates[,c("rate4site", "jalview", "site")]
plot.cons$rate4site <- plot.cons$rate4site*-1 + max(plot.cons$rate4site) #best to transform the rate4site measure. not very intuative
plot.cons.long <- melt(plot.cons, id.vars=c("site"))


tmp <- plot.cons[plot.cons$rate4site > 3,]
range <- data.frame(ranges = findIntRuns(tmp$site), stringsAsFactors = F)
range$xmin <- str_split(range$ranges, "-", simplify = T)[,1]
range$xmax <- str_split(range$ranges, "-", simplify = T)[,2]
range[range$xmax == "",]$xmax <- range[range$xmax == "",]$xmin
range$xmin <- as.numeric(range$xmin)
range$xmax <- as.numeric(range$xmax)

p.rate4site <- ggplot() +
  geom_rect(data=range, mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), alpha=0.5) +
  geom_line(data=plot.cons, aes(x=site, y=rate4site), color="grey40", size=1.1) +
  geom_ribbon(data=plot.cons, aes(ymin=rate4site,x=site), ymax=Inf, fill="white") +
  theme_classic() +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none") +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="#018571", alpha=0.75)
p.rate4site

tmp2 <- plot.cons[plot.cons$jalview > 6.6,]  ##cutoff about the same
range2 <- data.frame(ranges = findIntRuns(tmp2$site), stringsAsFactors = F)
range2$xmin <- str_split(range2$ranges, "-", simplify = T)[,1]
range2$xmax <- str_split(range2$ranges, "-", simplify = T)[,2]
range2[range2$xmax == "",]$xmax <- range2[range2$xmax == "",]$xmin
range2$xmin <- as.numeric(range2$xmin)
range2$xmax <- as.numeric(range2$xmax)

p.jalview <- ggplot() +
  geom_rect(data=range2, mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), alpha=0.5) +
  geom_line(data=plot.cons, aes(x=site, y=jalview), color="grey40", size=1.1) +
  geom_ribbon(data=plot.cons, aes(ymin=jalview,x=site), ymax=Inf, fill="white") +
  theme_classic() +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none") +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="#018571", alpha=0.75) +
  ylab("jalview")
p.jalview

foldx2$site <- as.numeric(foldx2$site)
p.fold <- ggplot(foldx2, aes(x=site, y=fold)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="#018571", alpha=0.75) +
  geom_hline(yintercept = 0, color="grey") +
  geom_jitter(height = 0.1, width=0.3, alpha=0.3) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ylab(expression(Delta~Delta~G[fold])) +
  xlab("") +
  theme(text = element_text(size=15), axis.title.x = element_blank(), axis.text.x = element_blank())
p.fold

p.bind <- ggplot(foldx2, aes(x=site, y=trimer)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="#018571", alpha=0.75) +
  geom_hline(yintercept = 0, color="grey") +
  geom_jitter(height = 0.1, width=0.3, alpha=0.3) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ylab(expression(Delta~Delta~G[trimer])) +
  xlab("") +
  theme(text = element_text(size=15), axis.title.x = element_blank())
p.bind


plot_grid(p.rate4site, p.jalview, p.fold, p.bind, nrow=4, align = 'v')
ggsave(filename = "plots/combo.pdf", width = 8.5, height= 8.5, units = "in")

#secondary structure from uniprot
#apparently only beta strands are known
sec_str <- data.frame(sites=c(13:17, 26:28, 30:33, 36:49, 51:61, 69:83, 88:98, 108:110, 115:117, 120:130, 134:137, 139:150, 152:163),
                      y=1, lab="B-sheets")
sec_str2 <- data.frame(xmin=c(13,26,30,36,51,69,88,108,115,120,134,139,152),
                       xmax=c(17,28,33,49,61,83,98,110,117,130,137,150,163))
sec_str2$ymin <- 0
sec_str2$ymax <- 1


###found interaction stuff in McKenna et al. 1994. J. Mol. Biol. (1994) 237, 517-543
##ggint from table 2 (2-fold interactions) and table 5 (5-fold interactions) 
gg_int2 <- data.frame(sites=c(1:8, 10, 13, 15, 70, 78, 84, 86:87, 89, 108:111, 113:117, 120, 122, 124:125, 134, 148, 156, 167:167, 171),
                      y="2", lab="Protein interfaces")
gg_int <- data.frame(xmin=c(1,10,13,15,70,78,84,86,89,108,113,120,122,124,134,148,156,167,171),
                       xmax=c(8,10,13,15,70,78,84,87,89,111,117,120,122,125,134,148,156,168,171))
gg_int$ymin <- 1
gg_int$ymax <- 2

ggplot(rsa_sum, aes(x=RESIDUE, y=RSA2)) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 0.1) +
  geom_line()+
  theme_classic()

#surface <- data.frame(rsa_high=rsa_sum[rsa_sum$RSA2 > 0.5, ]$RESIDUE, rsa_high_stop=rsa_sum[rsa_sum$RSA2 > 0.5, ]$RESIDUE)
surface <- data.frame(sites=rsa_sum[rsa_sum$RSA2 > 0.3, ]$RESIDUE, y=3, lab="Surface")

buried <- data.frame(sites=rsa_sum[rsa_sum$RSA2 < 0.05, ]$RESIDUE, y=4, lab="Buried")
#names(surface) <- c("xmin","xmax")
#surface$ymin <- 2
#surface$ymax <- 3

#p.str <- ggplot() +
#  geom_rect(data=sec_str2, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="#a6cee3", fill="#a6cee3") +
#  geom_rect(data=gg_int, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="#b2df8a", fill="#b2df8a") +
#  geom_rect(data=surface, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="#c2a5cf", fill="#c2a5cf") +
#  theme_classic() +
#  xlab("") +
#  theme(text = element_text(size=15), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
#  xlim(c(0,175))
#p.str

p.str <- ggplot() +
  geom_tile(data=sec_str, mapping=aes(x=sites, y=lab), color="#a6cee3", fill="#a6cee3") +
  geom_tile(data=gg_int2, mapping=aes(x=sites, y=lab), color="#b2df8a", fill="#b2df8a") +
  geom_tile(data=surface, mapping=aes(x=sites, y=lab), color="#c2a5cf", fill="#c2a5cf") +
  geom_tile(data=buried, mapping=aes(x=sites, y=lab), color="#7b3294", fill="#7b3294") +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="#018571", alpha=0.75) +
  theme_classic() +
  xlab("") +
  theme(text = element_text(size=15), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  xlim(c(0,175)) +
  scale_y_discrete(limits=c(  "B-sheets", "Protein interfaces","Surface", "Buried"))
p.str

plot_grid(p.rate4site, 
          p.jalview, 
          p.fold, 
          p.bind + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
          p.str, nrow=5, align = 'v')
ggsave(filename = "plots/combo_dists.pdf", width = 10, height= 9, units = "in")





####################################################################################################################################################################################
#plot differences in plaque size between different codons
datCodons <- read.csv("results/data_combined_allcods.csv", stringsAsFactors = F)
 codDiffs <- 
  datCodons %>%
  #filter(!codon == "XXX") %>%
  filter(!is.na(start.cod)) %>%
  group_by(site, aa) %>%
  filter(n() > 1)

##renormalize with dropping unobserved codons. plot looks funny b/c lots of expected values that are not observed
codDiffs[is.na(codDiffs$aacount),]$codonfreq <- NA
codDiffs$start.cod.obs <- codDiffs$start.cod
codDiffs[which(is.na(codDiffs$codonfreq)),]$start.cod.obs <- 0
codDiffs[which(is.na(codDiffs$codonfreq)),]$codonfreq <- 0

codDiffs <- codDiffs %>%
  group_by(site) %>%
  mutate(norm = start.cod.obs/sum(start.cod.obs))
codDiffs$start.cod.obs <- codDiffs$norm
codDiffs <- codDiffs[,-which(names(codDiffs)=="norm")]

codDiffs$codonOE <- codDiffs$codonfreq-codDiffs$start.cod
codDiffs$codonOE <- (codDiffs$codonfreq-codDiffs$start.cod.obs)/ codDiffs$start.cod.obs ##values calculated only for observed

#get wt codons
cods <- data.frame(site = 1:176, seq = substring(ref, seq(1,nchar(ref),3), seq(3,nchar(ref),3)))
codDiffs$wt <- cods[codDiffs$site,]$seq
codDiffs$wt <- codDiffs$wt == codDiffs$codon
##issue with one sample that only had 1 codon show up
codDiffs[codDiffs$sample.names == "4-161",]$codonOE <- 0
codDiffs[codDiffs$sample.names == "4-54",]$codonOE <- 0


limit <- max(abs(na.omit(codDiffs$codonOE))) * c(-1,1)
ggplot(codDiffs,aes(x=as.factor(site), y=codon, fill=codonOE, color=wt)) +
  geom_tile(width=0.95, height=0.95, size=0.5) +
  theme_classic() +
  #scale_fill_distiller("O-E/E", palette = "RdYlBu",na.value = "white") +
  scale_fill_distiller("O-E/E", palette = "RdYlBu", limit=limit, na.value = "white") +
  scale_color_manual(values=c("white", "black"), guide="none") +
  #scale_fill_distiller("Observ-Expect", palette = "RdYlBu", na.value = "white") +
  #geom_tile(data=cods.short, aes(x=as.factor(site), y=seq), fill="black")+
  facet_grid(aa~.,scales = "free_y", space="free_y", switch="y", drop = T) +
  theme(strip.placement = "outside",                
        strip.background = element_rect(fill = "white", colour = "white"), 
        axis.title = element_blank())
ggsave(filename = "plots/codon_freq.pdf", width = 6.5, height= 9.5, units = "in")


# 
# 
# 
# #####plot distance matrix heatmap
# library(distances)
# library(gplots)
# dat <- read.table(file="results/fg_alpha", stringsAsFactors = F)
# dat.s <- data.frame(dat[,7:9])
# names(dat.s) <- c("x","y","z")
# #dist.mat <- matrix(nrow=nrow(dat),ncol=nrow(dat))
# 
# euclidean_distance <- function(p,q){
#   sqrt(sum((p - q)^2))
# }
# 
# euclidean_distance(dat.s[1,],dat.s[2,])
# 
# #dist.mat <- stats::dist(dat.s[1:10,], diag = T, upper=T) ###huh. having a problem with using all the data
# 
# #trying differenc package
# dist.mat.2 <- distances(dat.s)
# dist.mat.2 <- as.matrix(distances(dat.s))
# ggplot(dist.mat.2) +
#   geom_tile
# 
# heatmap.2(dist.mat.2, Rowv=F, Colv=F, dendrogram = "none", labRow=F, labCol = F, trace="none", col=c("blue", "lightgrey", "white"))
# 
# 
# colors = c(seq(0,10,length=100),seq(10.01,20,length=100),seq(20.01,105,length=100))
# my_palette <- colorRampPalette(c("red", "black", "skyblue"))(n = 299)
# heatmap.2(dist.mat.2, col=my_palette, 
#           breaks=colors, density.info="none", trace="none", Rowv=F, Colv=F, 
#           dendrogram="none", symm=F,symkey=F,symbreaks=T, scale="none",
#           labRow=F, labCol=F)
# 
# dist.mat.g <- dist.mat.2[427:601,427:601]
# colnames(dist.mat.g) <- 1:175
# row.names(dist.mat.g) <- 1:175
# dist.mat.g.long <- melt(dist.mat.g)
# dist.mat.g.long$lin_dist <- abs(dist.mat.g.long$Var1-dist.mat.g.long$Var2)
# long.dist.g <- dist.mat.g.long[dist.mat.g.long$lin_dist > 5,]
# 
# p.dist <- ggplot(long.dist.g, aes(x=Var1, y=value)) +
#   geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
#   #geom_hline(yintercept = 0, color="grey") +
#   geom_jitter(height = 0.1, width=0.3, alpha=0.3) +
#   stat_summary(fun.y=median, geom="point", shape=18,
#                size=2, color="red") +
#   theme_classic() +
#   theme(axis.text.x=element_text(angle=90,hjust=1)) +
#   ylab(expression(Delta~Delta~G[trimer])) +
#   xlab("") +
#   theme(text = element_text(size=15), axis.title.x = element_blank())
# p.bind
# 
# 
# tmp <- long.dist.g[,c("Var1", "value")]
# library(dplyr)
# dist.agg <- tmp %>% group_by(Var1) %>%
#   summarise_each(funs(.[which.max(abs(.))]))
# 
# 
# p.dist.agg <- ggplot(dist.agg, aes(x=Var1, y=value)) +
#   geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
#   geom_line(color="grey40", size=1.1) +
#   theme_classic() +
#   ylab("avg atomic dist") +
#   theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none")
# p.dist.agg
# 
# plot_grid(p.rate4site, 
#           p.jalview, 
#           p.fold, 
#           p.bind + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
#           p.dist.agg,
#           p.str, nrow=6, align = 'v')
# ggsave(filename = "plots/combo_dists.pdf", width = 10, height= 10, units = "in")
