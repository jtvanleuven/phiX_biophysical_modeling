calc.posteriors.details <- function(data, wt, viable.01=NA, cat.probs.here){
  post.matrix <- as.data.frame(matrix(nrow=1, ncol=9))
  colnames(post.matrix) <- c("N", "n.true","n.obs", "post.true", "post.obs", "post.mode", "mode.prob", "up.95", "up.99")
  
  n.true <- NA
  n <- sum(data)
  zeros <- which(data==0)
  if (length(zeros>0)){
    obs.freqs <- data/sum(data)
    hyp <- as.numeric(data>0)
    n.obs <- sum(hyp)
    n.unobs <- length(cats)- n.obs
    obs.AAs <- cats[which(hyp==1)]
    unobs.AAs <- cats[which(hyp==0)]
    hyp.01 <- paste(hyp, collapse="")
    
    
    #generate alternative hypotheses
    hyp.vects <- vector("list", n.unobs)
    names(hyp.vects) <- unobs.AAs
    for (i in 1:n.unobs){
      hyp.vects[[i]] <- c(0,1)
    }
    alt.hyps <- expand.grid(hyp.vects)
    
    viable.hyps <- matrix(nrow=length(alt.hyps[,1]), ncol=length(cats))
    colnames(viable.hyps) <- cats
    for (i in 1:length(obs.AAs)){
      col <- which(cats==obs.AAs[i])
      viable.hyps[,col] <- 1
    }
    for (i in 1:length(unobs.AAs)){
      col <- which(cats==unobs.AAs[i])
      viable.hyps[,col] <- alt.hyps[,i]
    }
    
    if (data[wt]==0){  # If wildtype (which must be viable) is unobserved 
      throwouts <- which(viable.hyps[,wt]==0)
      vaible.hyps <- viable.hyps[-throwouts,]
    }
    viab.hyps.01 <- apply(viable.hyps, 1, function(x) paste(x, collapse=""))
    obs.hyp <- which(viab.hyps.01==hyp.01)
    if (is.na(viable.01)==FALSE){
      true.hyp <- which(viab.hyps.01==viable.01)
      n.true <- sum(as.numeric(unlist(strsplit(viable.01, split=""))))
    } else{
      true.hyp <- NA
      n.true <- NA
    }
    
    probs.by.via.hyp <- t(apply(viable.hyps, 1, function(x) x*cat.probs.here))
    probs.by.via.hyp <- t(apply(probs.by.via.hyp, 1, function(x) x/sum(x)))
    p.data.by.via.hyp <- t(t(apply(probs.by.via.hyp, 1, function(x) dmultinom(data, prob=x))))
    
    #per.hyp.prob <- 1/(2^length(cats))
    #bayes.num <- p.data.hyp*per.hyp.prob
    #bayes.denom <- sum(p.data.by.via.hyp*per.hyp.prob)
    #posterior <- bayes.num/bayes.denom
    
    posteriors <- sapply(p.data.by.via.hyp, function(x) x/sum(p.data.by.via.hyp))
    post.obs.prob <- posteriors[obs.hyp]
    if (is.na(true.hyp)){
      post.true.prob <- NA
    } else{
      post.true.prob <- posteriors[true.hyp]  
    }
    post.by.nviab <- rep(0, 20)
    n.by.viab.hyps <- apply(viable.hyps, 1, function(x) sum(x))
    post.by.nviab <- sapply(seq(1,20), function(x) sum(posteriors[which(n.by.viab.hyps==x)]))
    post.mode <- seq(1,20)[which.max(post.by.nviab)]
    post.mode.prob <- post.by.nviab[post.mode]
    up.bound.95 <- seq(1,20)[min(which(cumsum(post.by.nviab)>0.95))]
    up.bound.99 <- seq(1,20)[min(which(cumsum(post.by.nviab)>0.99))]
    
    # --- Posteriors by AA ---
    rows.by.AA <- lapply(seq(1,20), function(x) which(viable.hyps[,x]==1))
    prob.by.AA <- unlist(lapply(rows.by.AA, function(x) sum(posteriors[x])))
    names(prob.by.AA) <- colnames(viable.hyps)
    
    
  } else{
    post.obs.prob <- 1
    post.true.prob <- 1
    post.by.nviab <- c(rep(0,19),1)
    n.obs <- 20
    post.mode <- n.obs
    post.mode.prob <- 1
    up.bound.95 <- n.obs
    up.bound.99 <- n.obs
    prob.by.AA <- rep(1, 20)
    names(prob.by.AA) <- cats
  }
  j <- 1
  post.matrix$N[j] <- n
  post.matrix$n.true[j] <- n.true
  post.matrix$n.obs[j] <- n.obs
  post.matrix$post.true[j] <- post.true.prob
  post.matrix$post.obs[j] <- post.obs.prob
  post.matrix$post.mode[j] <- post.mode
  post.matrix$mode.prob[j] <- post.mode.prob
  post.matrix$up.95[j] <- up.bound.95
  post.matrix$up.99[j] <- up.bound.99
  
  return(list(post.matrix=post.matrix,  post.by.nviab=post.by.nviab, prob.by.AA=prob.by.AA))
}