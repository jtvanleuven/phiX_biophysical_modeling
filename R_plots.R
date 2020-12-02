library(ggplot2)
library(stringr)
library(cowplot)
library(Biostrings)
library(seqinr)
library(reshape2)

ref <- readDNAStringSet("rawdata/G_wt.fasta")
refaa <- translate(ref)

##problems
###depends on plaque size code from old stuff still. currently reads in "plot.mat.csv" which is created in "R_visualize_data_v2.R"


plot.mat <- read.csv("results/plot.mat.csv", row.names = 1, stringsAsFactors = F)  ###need to add code to the repository that generates this file
foldx <- read.csv(file="results/ddG_foldX_GProtein.csv", header=T, stringsAsFactors = F)
foldx2 <- foldx
foldx2$aa <- str_sub(foldx2$mut,1,1)
foldx2$var <- str_sub(foldx2$mut,-1,-1)
foldx2$site <- str_extract(foldx2$mut,'\\d+')

##pretty sure craig made this table
dat <- read.csv("results/Combined_pyR_FoldX_data_v3.csv", row.names = 1, stringsAsFactors = F)
##pretty sure I made this table. Contains plaque size and codon recovery data.
datCodons <- read.csv("results/data_combined.csv", stringsAsFactors = F)

seq.counts <- read.csv(file="results/protein_G_mutagenesis_AA_counts.csv", header=T)
##F-G data
GFbind <- read.csv("results/G_mut_GF_bind_ddg_stdev.csv", header = F, stringsAsFactors = F)

dat <- read.csv("results/Combined_pyR_FoldX_data_v3.csv", row.names = 1, stringsAsFactors = F)


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

p.rate4site <- ggplot(plot.cons, aes(x=site, y=rate4site)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
  geom_line(color="grey40", size=1.1) +
  theme_classic() +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none")
p.rate4site

p.jalview <- ggplot(plot.cons, aes(x=site, y=jalview)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
  geom_line(color="grey40", size=1.1) +
  theme_classic() +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none") +
  ylab("jalview")
p.jalview

foldx2$site <- as.numeric(foldx2$site)
p.fold <- ggplot(foldx2, aes(x=site, y=fold)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
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
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
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
sec_str <- c(13:17, 26:28, 30:33, 36:49, 51:61, 69:83, 88:98, 108:110, 115:117, 120:130, 134:137, 139:150, 152:163)
sec_str2 <- data.frame(xmin=c(13,26,30,36,51,69,88,108,115,120,134,139,152),
                       xmax=c(17,28,33,49,61,83,98,110,117,130,137,150,163))
sec_str2$ymin <- 0
sec_str2$ymax <- 1


###found interaction stuff in McKenna et al. 1994. J. Mol. Biol. (1994) 237, 517-543
##ggint from table 2 (2-fold interactions) and table 5 (5-fold interactions) 
gg_int <- data.frame(xmin=c(1,10,13,15,70,78,84,86,89,108,113,120,122,124,134,148,156,167,171),
                       xmax=c(8,10,13,15,70,78,84,87,89,111,117,120,122,125,134,148,156,168,171))
gg_int$ymin <- 1
gg_int$ymax <- 2

p.str <- ggplot() +
  geom_rect(data=sec_str2, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="#a6cee3", fill="#a6cee3") +
  geom_rect(data=gg_int, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="#b2df8a", fill="#b2df8a") +
  theme_classic() +
  xlab("") +
  theme(text = element_text(size=15), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlim(c(0,175))
p.str

plot_grid(p.rate4site, 
          p.jalview, 
          p.fold, 
          p.bind + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
          p.str, nrow=5, align = 'v')


#####plot distance matrix heatmap
library(distances)
library(gplots)
dat <- read.table(file="results/fg_alpha", stringsAsFactors = F)
dat.s <- data.frame(dat[,7:9])
names(dat.s) <- c("x","y","z")
#dist.mat <- matrix(nrow=nrow(dat),ncol=nrow(dat))

euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}

euclidean_distance(dat.s[1,],dat.s[2,])

#dist.mat <- stats::dist(dat.s[1:10,], diag = T, upper=T) ###huh. having a problem with using all the data

#trying differenc package
dist.mat.2 <- distances(dat.s)
dist.mat.2 <- as.matrix(distances(dat.s))
ggplot(dist.mat.2) +
  geom_tile

heatmap.2(dist.mat.2, Rowv=F, Colv=F, dendrogram = "none", labRow=F, labCol = F, trace="none", col=c("blue", "lightgrey", "white"))


colors = c(seq(0,10,length=100),seq(10.01,20,length=100),seq(20.01,105,length=100))
my_palette <- colorRampPalette(c("red", "black", "skyblue"))(n = 299)
heatmap.2(dist.mat.2, col=my_palette, 
          breaks=colors, density.info="none", trace="none", Rowv=F, Colv=F, 
          dendrogram="none", symm=F,symkey=F,symbreaks=T, scale="none",
          labRow=F, labCol=F)

dist.mat.g <- dist.mat.2[427:601,427:601]
colnames(dist.mat.g) <- 1:175
row.names(dist.mat.g) <- 1:175
dist.mat.g.long <- melt(dist.mat.g)
dist.mat.g.long$lin_dist <- abs(dist.mat.g.long$Var1-dist.mat.g.long$Var2)
long.dist.g <- dist.mat.g.long[dist.mat.g.long$lin_dist > 5,]

p.dist <- ggplot(long.dist.g, aes(x=Var1, y=value)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
  #geom_hline(yintercept = 0, color="grey") +
  geom_jitter(height = 0.1, width=0.3, alpha=0.3) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ylab(expression(Delta~Delta~G[trimer])) +
  xlab("") +
  theme(text = element_text(size=15), axis.title.x = element_blank())
p.bind


tmp <- long.dist.g[,c("Var1", "value")]
library(dplyr)
dist.agg <- tmp %>% group_by(Var1) %>%
  summarise_each(funs(.[which.max(abs(.))]))


p.dist.agg <- ggplot(dist.agg, aes(x=Var1, y=value)) +
  geom_vline(xintercept = as.numeric(unique(datCodons$site)), color="grey", alpha=0.5) +
  geom_line(color="grey40", size=1.1) +
  theme_classic() +
  ylab("avg atomic dist") +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none")
p.dist.agg

plot_grid(p.rate4site, 
          p.jalview, 
          p.fold, 
          p.bind + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
          p.dist.agg,
          p.str, nrow=6, align = 'v')
ggsave(filename = "plots/combo_dists.pdf", width = 10, height= 10, units = "in")
