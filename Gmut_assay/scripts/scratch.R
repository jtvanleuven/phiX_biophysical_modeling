sample_ls <- list()
n_singletons <- list()
mismatches <- list()
for (s in 1:length(samples)){
  if (nrow(zoelist[[s]]$mut_info)==1){
    sample_ls[[s]] <- samples[s]
    n_singletons[[s]] <- nrow(zoelist[[s]]$singletons)
    mismatches[[s]] <- zoelist[[s]]$singletons$mismatches
  }
}


vars <- unique(yesol_df2$AAsub)
vars_df <- as.data.frame(vars)
samples <- unique(yesol_df2$sample)
for (i in 1:length(samples)){
  df_i <- filter(yesol_df2, sample == samples[i])
  for (j in 1:length(vars)){
    if (vars[j] %in% df_i$AAsub){
      vars_df[j, samples[i]] <- df_i[df_i$AAsub==vars[j],]$count
    } else {
      vars_df[j, samples[i]] <- 0
    }
  }
}


sam_path <- "~/Projects/phix_host/PrelimAssay/data/dada2/final_sam/t35_rep1.sam"
z <- read.table(sam_path, skip = 2, stringsAsFactors = F)


offs <- yesol_df %>% filter(in_target == 0) 
yesol_df %>% 
  filter(in_target == 0) %>%
  ggplot(aes(x = site)) +
  geom_bar()

sampath <- "~/Projects/phix_host/PrelimAssay/data/dada2/final_sam/t35_rep1.sam"
sam <- read.table(sampath, skip = 2, stringsAsFactors = F)
ct <- str_split(sam[,1], "=", simplify = TRUE)[,2]
sam$seqID <- str_split(sam[,1], ";", simplify = T)[,1] #sq1, sq2, etc
sam$count <- as.numeric(str_replace(ct, ";", ""))
sam$mismatches <- as.numeric((str_split(sam[,12], ":", simplify = TRUE)[,3]))
# Save singletons (counts == 1) as separate df
singletons <- filter(sam, count==1) %>% select(V1, mismatches)
# delete singleton rows from main df
sam <- filter(sam, count > 1)
# Loop through each read, compare codons
# first, get ref codons
wt.seq <- sam[sam$mismatches==0,]$V10
start <- seq(1, nchar(as.character(wt.seq)), 3)
stop <- start+2
wt.cod <- str_sub(wt.seq, start, stop)

par(mfrow=c(1,2))
hist(yesol_df$site)
hist(zoe_df$site)


df1 <- select(yesol_avgscores, vars, avg_score)
df2 <- select(zoe_avgscores, vars_z, avg_score)
dfjoin <- inner_join(df1, df2, by = c("vars" = "vars_z"))


### Figure out synonymous muts
# get glib and wt.cod from prev scripts (use path below)

sampath <- "data/sed_sam/t0.dada.sam"
sam <- read.table(sampath, stringsAsFactors = F)

# make a df listing all AA subs in the primers list that are synonymous muts
sm <- select(glib, mut.corr, pos, cod.seq) %>%
  mutate(AAfrom = sapply(mut.corr, function(x) str_split(gsub("\\d","",x),"", simplify = T)[1]),
         AAto = sapply(mut.corr, function(x) str_split(gsub("\\d","",x),"", simplify = T)[2])) %>%
  filter(AAfrom == AAto)

# get wt codons to compare
wt <- as_tibble(wt.cod) %>%
  mutate(pos = seq(1:176))
colnames(wt) <- c("wt.cod", "pos")

# match by position
matched <- right_join(sm, wt, by = "pos") %>%
  select(mut.corr, pos, cod.seq, wt.cod) 
colnames(matched) <- c("muts", "pos", "mut.cod", "wt.cod")

write.csv(matched, "~/Desktop/synonymous_muts.csv", row.names = F)



# 20200723
viability <- read.csv("~/Projects/phix_host/phix_jt/Craigs_21Site_ddG_Viability_data.csv")[,-c(1:2)]
viability <- select(viability, mutation, obs)
colnames(viability) <- c("vars","obs")
joined_df <- left_join(combined_scores, viability, by = "vars")
joined2 <- left_join(viability, combined_scores, by = "vars")
joined3 <- inner_join(viability, combined_scores, by = "vars")

# regress binary viability on fitness score
reg <- glm(obs~betaML, data = joined3, family = binomial)

# plot
x <- seq(-0.11, 0.07, 0.01)
y <- predict(reg, list(betaML = x), type="response")
plot(joined3$obs~joined3$betaML)
lines(x, y)


# 20200728
yesol_vars_wd[rowSums(yesol_vars_wd[ ,5:17]) == 0, "vars"]






combined_scores[combined_scores$vars == "H58S",]
yesol_vars_wd %>%
  filter(vars == "H58S") %>% select(t0:t35_rep3, t70_rep1:t70_rep3)
yesol_vars_wd %>%
  filter(vars == "T105Q,V164A") %>% select(t0:t35_rep3, t70_rep1:t70_rep3)
combined_scores[combined_scores$vars == "T105Q,V164A",]

