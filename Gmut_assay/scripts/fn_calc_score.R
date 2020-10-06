
# Input: vars df, sample times as a vector, reps as numeric
#         df should be in a long form with cols: vars, count, t(0, 30,..), rep("rep1","rep2",...)
# output: df with vars, avg_score, sd, site, AA from and to
calc_score <- function(df, t, reps){
  vars <- unique(df$vars)
  score_df <- as.data.frame(matrix(data = NA, nrow = length(vars), ncol = 1+reps))
  se_df <- as.data.frame(matrix(data = NA, nrow = length(vars), ncol = 1+reps))
  pval_df <- as.data.frame(matrix(data = NA, nrow = length(vars), ncol = 1+reps))
  score_df[,1] <- vars
  se_df[,1] <- vars
  pval_df[,1] <- vars
  # get numbers for t0 (same for all reps, so doesn't need to be in the loop)
  t0 <- filter(df, time == 0)
  Cwt0 <- t0[t0$vars == "wt",]$count
  Mv0 <- sapply(t0$count, function(x) log((x+.5)/(Cwt0+.5)))
  Vv0_inv <- sapply(t0$count, function(x) 1/(1/(x+.5)+(1/(Cwt0+0.5))))
  rep_names <- vector() # to name columns later
  for (i in 1:reps){
    rep_i <- paste("rep", i, sep = "")
    rep_names <- c(rep_names, rep_i)
    dfs <- list()
    Cwts <- list()
    Mvs <- list()
    Vv_invs <- list()
    dfs[[1]] <- t0
    Cwts[[1]] <- Cwt0
    Mvs[[1]] <- Mv0
    Vv_invs[[1]] <- Vv0_inv
    for (j in 2:length(t)){
      dfs[[j]] <- filter(df, time == t[j], rep == rep_i)
      Cwts[[j]] <- dfs[[j]][dfs[[j]]$vars == "wt", ]$count
      Mvs[[j]] <- sapply(dfs[[j]]$count, function(x) log((x+.5)/(Cwts[[j]]+.5)))
      Vv_invs[[j]] <- sapply(dfs[[j]]$count, function(x) 1/(1/(x+.5)+(1/(Cwts[[j]]+0.5))))
    }
    for (k in 1:nrow(score_df)){
      Ms <- vector()
      wts <- vector()
      for (l in 1:length(t)){
        Ms <- c(Ms, Mvs[[l]][k])
        wts <- c(wts, Vv_invs[[l]][k])
      }
      reg <- lm(Ms~t, weights = wts)
      score_df[k, i+1] <- reg$coefficients["t"]
      se_df[k, i+1] <- summary(reg)$coefficients["t", "Std. Error"]
      pval_df[k, i+1] <- summary(reg)$coefficients["t", "Pr(>|t|)"]
    }
  }
  colnames(score_df) <- c("vars", rep_names)
  colnames(se_df) <- c("vars", rep_names)
  colnames(pval_df) <- c("vars", rep_names)
  out <- list()
  out[["scores"]] <- score_df
  out[["se"]] <- se_df
  out[["pval"]] <- pval_df
  return(out)
}
