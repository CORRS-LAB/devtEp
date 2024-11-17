# TODO

1. This project aims to identify cell type specific DNA methylation regions associated with phenotypes using RRBS data.
   To derive a robust estimation of immune cell fractions, we first perform a benchmarking study to compare 19 deconvolution algorithms using RRBS data.
   Four key steps are included in this project: marker selection, benchmarking of deconvolution algorithms, estimating the immune cell fractions from real RRBS data, identifing cell type specific DNA methylation regions associated with phenotypes.

2. Immune Cell Type Specific DNA Methylation Regions associate with phenotypes
   
3. HELP Manual:
    - marker selection and preparing the reference signatures
    - benchmarking of deconvolution algorithms
    - estimating the immune cell fractions from real RRBS data
    - identifying cell type specific DNA methylation regions associated with phenotypes


# Demo: marker selection and preparing the reference signatures.

Regions were selected as deconvolution makers using the algorithm described by Luo et al.[PMID: 22110609; PMID: 23284283]. 
To increase the computational efficiency, we optimized the original marker selection algorithm using C/C++ language and incorporated it into our R package “devtEp”. 
For each region, two-sample t-tests were applied to assess methylation difference between the target cell type and all other cell types. 
The top 200 regions with lowest P values were selected for each cell type, and from these, the 100 regions with the highest mean methylation differences were retained. 
After removing duplicated regions across cell types, the final reference signature was derived.




```{r}
library(devtEp)
library(plyr)

## A function to select top N values in a data matrix, return the row ids.
select_top_N <- function(dt.input,N, dec){

   # Find the top N values in each column
   top_N_indices <- apply(dt.input, 2, function(x) order(x, decreasing = dec)[1:N])

   # Select the rows with the top N values in each column
   selected_rows <- data.matrix(unique(unlist(top_N_indices)))
   return(selected_rows)
}


data <- demo_devEpt_1k

names <- colnames(data)  
colnames(data)[-c(1:3)] <- paste(names[-c(1:3)],".Z00R", sep = "") #### TO DO: Can we remove this?
uu <- selectSite(df = data, skip.col = 3)

n <- 100
row.ids <- c(unlist(select_top_N(uu$pvalue,2*n, FALSE)))
row.ids <- row.ids[order(row.ids)]

Dif <- uu$abs_diff[row.ids,]
rownames(Dif) <- row.ids

top.dif <- apply(Dif, 2, function(x) order(x, decreasing = TRUE)[1:n])
ids <- unique(as.numeric(rownames(Dif[top.dif,])))
row.out <- ids[order(ids)]

out <- data[row.out,]
colnames(out) <- names
write.table(out, "reference.mker.ref",row.names = F, sep = "\t", quote = F)

```


# Demo: benchmarking of deconvolution algorithms

```{r}
### Functions used
library(reshape)
library(glmnet)
library(bvls)
library(EMeth)
library(nnls)
library(preprocessCore)
library(quadprog)
library(stringr)
library(EpiDISH)
library(minfi)
library(EMeth)
library(Matrix)
library(MethylResolver)
library(TOAST)
library(corpcor)
library(writexl)
library(limma)
library(ADAPTS)
library(ICeDT)
library(FARDEEP)
library(EPIC)

####minfi
projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = FALSE) {
  if (is.null(contrastCellType)) {
    Xmat <- coefCellType
  } else {
    Xmat <- tcrossprod(coefCellType, contrastCellType)
  }
  
  nCol <- dim(Xmat)[2]
  if (nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(
      apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  } else {
    nSubj <- dim(Y)[2]
    
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    
    if (nonnegative) {
      if (lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      } else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve.QP(
          Dmat = Dmat,
          dvec = crossprod(Xmat[obs,], Y[obs,i]),
          Amat = Amat,
          bvec = b0vec)$sol
      }
    } else {
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    mixCoef
  }
}

## This is the normalization fucntion; dataset is in data.matrix format.
normalization <- function(u, dataset){
  if(u == 'none'){
    #None
    method <- 'none'
    normaliz <- function(x){
      return(x)
    }
    
  }else if (u == 'z_score'){
    #Z-score
    method = 'z_score'
    normaliz <- function(x){
      return((x-mean(na.omit(x)))/sd(na.omit(x)))
    }
    
  } else if (u == 'min_max'){
    #Min-max
    method <- 'min_max'
    normaliz <- function(x){
      return((x-min(na.omit(x)))/(max(na.omit(x))-min(na.omit(x))))
    }
    
  } else if (u == 'col_z_score'){
    #Column Z-score
    method <- 'col_z_score'
    normaliz <- function(x){
      cols <-function(y){
        return((y-mean(na.omit(y)))/sd(na.omit(y)))
        
      }
      apply(x, 2, cols)
      
    }
    
  } else if (u == 'col_min_max'){
    #Column min-max
    method <- 'col_min_max'
    
    normaliz <- function(x){
      cols <-function(y){
        return((y-min(na.omit(y)))/(max(na.omit(y)-min(na.omit(y)))))
        
      }
      apply(x, 2, cols)
      
    }
    
  } else if (u == 'QN'){
    #Quantile normalization
    
    method <- 'QN'
    normaliz <- function(x){
      x <- na.omit(x)
      rows <- rownames(x)
      cols <- colnames(x)
      x <- normalize.quantiles(x)
      rownames(x) <- rows
      colnames(x) <- cols
      return(x)
    }
    #
  }  else if (u == 'lognorm'){
    method <- 'lognorm'
    normaliz <- function(x){
      x <- log(na.omit(x)+1)
      x[abs(x)==Inf] <- NA
      return(na.omit(x))
    }
  }
  out_dataset <- normaliz(dataset)
  return(out_dataset)
}

prop_out <- function(prop, predicted_cell_proportions,method){
  
  colnames(predicted_cell_proportions) <- colnames(prop)
  predict.prop <- data.frame(sample = rownames(prop), predicted_cell_proportions)
  predict.out <- melt(predict.prop)
  colnames(predict.out)[3] <- "Predict proportion"
  
  real.prop <- data.frame(sample = rownames(prop), prop)
  real.out <- melt(real.prop)
  colnames(real.out)[3] <- "True proportion"
  
  tmp <- data.frame(merge(predict.out, real.out, by =c("sample", "variable")),method)
  return(tmp)
}

#### This is the deconvolution function.
deconvolution.compare <- function(df, ref, prop){
  out <- NULL
  #prop <- props
  real_cell_proportions <- prop
  
  #EMeth-Normal-----
  method = 'EMeth-Normal'
  maximum_nu = 0
  maximum_iter = 50
  predicted_cell_proportions <- matrix(NA, ncol = dim(prop)[2], nrow = dim(prop)[1])
  predicted_cell_proportions <- emeth(Y = as.matrix(df),eta = c(rep(0, dim(as.matrix(df))[2])), mu = as.matrix(ref), aber = FALSE, V = 'c',init = "default", family = "normal", nu = maximum_nu, maxiter = maximum_iter, verbose = TRUE)$rho
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #EMeth-Laplace-----
  method = 'EMeth-Laplace'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- emeth(Y = as.matrix(df),eta = c(rep(0, dim(as.matrix(df))[2])), mu = as.matrix(ref), aber = FALSE, V = 'c',init = "default", family = "laplace", nu = maximum_nu, maxiter = maximum_iter, verbose = TRUE)$rho
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #EMeth-Binomial-----
  method = 'EMeth-Binomial'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- emeth(Y = as.matrix(df),eta = c(rep(0, dim(as.matrix(df))[2])), mu = as.matrix(ref), aber = FALSE, V = 'b',init = "default", family = "normal", nu = maximum_nu, maxiter = maximum_iter, verbose = TRUE)$rho
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  
  #NNLS-----
  method = 'NNLS'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  for (i in 1:dim(real_cell_proportions)[1]){
    predicted_cell_proportions[i,] <- coef(nnls(as.matrix(ref), df[,i]))
  }
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method)) 
  
  #BLVS-----
  method = 'BLVS'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  for (i in 1:dim(real_cell_proportions)[1]){
    # predicted_cell_proportions[i,] <- coef(bvls(as.matrix(ref), df[rownames(ref),i], bl = c(rep(0, dim(ref)[2])), bu = c(rep(1, dim(ref)[2]))))
    predicted_cell_proportions[i,] <- coef(bvls(as.matrix(ref), df[,i], bl = c(rep(0, dim(ref)[2])), bu = c(rep(1, dim(ref)[2]))))
  }
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #Minfi-----
  method = 'Minfi'
  predicted_cell_proportions <- projectCellType(Y = as.matrix(df), coefCellType = as.matrix(ref))  
 
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #MethylResolver-----
  method = 'MethylResolver'
  methylMix <- as.data.frame(df)
  methylSig <- as.data.frame(ref)
  regressionFormula = as.formula(paste0("methylMix[,i] ~ ",paste(colnames(ref),sep="",collapse=" + ")))
  preds <- matrix(NA, ncol = length(colnames(ref)), nrow = length(colnames(df)))
  colnames(preds) <- colnames(ref)
  rownames(preds) <- colnames(df)
  j = 1
  for (i in colnames(df)){
    deconvoluteSample <- robustbase::ltsReg(regressionFormula, data = methylSig, alpha = 0.5)
    preds[j,] <- deconvoluteSample$coefficients[2:length(deconvoluteSample$coefficients)]
    j <- j+1
  }
  predicted_cell_proportions <- preds
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #OLS-----
  method = 'OLS'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  cmd = 'lm(formula = df[,i] ~ 0'
  for(i in colnames(ref)){
    cmd = paste0(cmd, paste0(' + ref[,"',i, '"]'))
  }
  cmd = paste0(cmd, ')')
  for (i in 1:dim(df)[2]){
    model <- eval(parse(text = cmd))
    m[i,] <- unlist(lapply(model$coefficients, function (x) x/sum(model$coefficients)))
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))

  #Ridge------
  method = 'Ridge'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  for (i in 1:dim(df)[2]){
    mod <- glmnet(x = as.matrix(ref),y=df[, i], alpha = 0)
    mod <- glmnet(x = as.matrix(ref),y=df[, i], alpha = 0, lambda = mod$lambda[which(mod$dev.ratio == max(mod$dev.ratio))])
    m[i,] <- coef(mod)[2:length(coef(mod))]
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #Elastic net------
  method = 'Elastic net'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  for (i in 1:dim(df)[2]){
    mod <- glmnet(x = as.matrix(ref),y=df[, i], alpha = 0.5)
    mod <- glmnet(x = as.matrix(ref),y=df[, i], alpha = 0.5, lambda = mod$lambda[which(mod$dev.ratio == max(mod$dev.ratio))])
    m[i,] <- coef(mod)[2:length(coef(mod))]
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #Lasso------
  method = 'Lasso'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  for (i in 1:dim(df)[2]){
    mod <- glmnet(x = as.matrix(ref),y=df[, i], alpha = 1)
    mod <- glmnet(x = as.matrix(ref),y=df[, i], alpha = 1, lambda = mod$lambda[which(mod$dev.ratio == max(mod$dev.ratio))])
    m[i,] <- coef(mod)[2:length(coef(mod))]
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #ICeDT------
  method = 'icedt'
  predicted_cell_proportions <- as.data.frame(t(ICeDT(as.matrix(df), as.matrix(ref), rhoConverge = 0.0015, maxIter_PP = 30)$rho))[,colnames(real_cell_proportions)]

  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))

  #FARDEEP-----
  method = 'FARDEEP'
  rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)))
  predicted_cell_proportions <- as.data.frame(fardeep(as.matrix(ref),as.matrix(df))$abs.beta)[,colnames(real_cell_proportions)]
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  
  #DCQ-----
  method = 'DCQ'
  predicted_cell_proportions <- as.data.frame(t(estCellPercent.DCQ(refExpr=as.matrix(ref), geneExpr=as.matrix(df))))[,colnames(real_cell_proportions)]/100
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #EpiDISH-CP------
  method = 'Epidish-CP'
  rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)))
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- epidish(as.matrix(df), as.matrix(ref),method= 'CP')$estF

  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))

  #Epidish-RPC-----
  method = 'Epidish-RPC'
  rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)))
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- epidish(as.matrix(df), as.matrix(ref),method= 'RPC')$estF

  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))

  #EpiDISH-CBS------
  method = 'Epidish_CBS'
  rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)))
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- epidish(as.matrix(df), as.matrix(ref),method= 'CBS')$estF

  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))

  #EPIC-----
  method = 'EPIC'
  rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)))
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  reflist <- list(
    refProfiles = as.matrix(ref),  # β值作为参考细胞类型的基因表达特征
    sigGenes = rownames(ref)
  )
  predicted_cell_proportions <- EPIC(bulk = as.matrix(df), reference = reflist)$cellFractions
  predicted_cell_proportions <- predicted_cell_proportions[,colnames(real_cell_proportions)]
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #TOAST-----
  method = 'TOAST'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- TOAST::myRefFreeCellMix(Y = as.matrix(df), mu0 =as.matrix(ref))$Omega
  predicted_cell_proportions <- predicted_cell_proportions[,colnames(real_cell_proportions)]
  
  out <- rbind(out, prop_out(prop, predicted_cell_proportions, method))
  
  #-----
  return(out)
}

####JSD计算函数
KL_divergence <- function(P, Q) {
  sum(P * log(P / Q))
}
calculate_jsd <- function(pred,true){
  epsilon <- 1e-10
  p <- ifelse(pred<0,0,pred)
  q <- ifelse(true<0,0,true)
  M <- 0.5 * (p + q)
  p_M <- KL_divergence(p + epsilon, M + epsilon)
  q_M <- KL_divergence(q + epsilon, M + epsilon)
  
  jsd <- 0.5 * (p_M + q_M)
  return(jsd)
}




region.proportion.estimation <- function(ref.file, bulk.file, prop.file, out.file, normalizations){
  
  result.out <- NULL
  
  for (u in normalizations){
    ref.dat = read.table(ref.file, header =T)
    ref.dt <- ref.dat[,-c(1:3)]
    ref.dt <- normalization(u, as.matrix(ref.dt))
    celltypes <- colnames(ref.dt) <- read.table(text =colnames(ref.dt), sep = ".")[,2]
    
    ref <- matrix(NA, ncol = length(unique(celltypes)), nrow = nrow(ref.dt))
    colnames(ref) <- unique(celltypes)
    for(i in colnames(ref)){
      ref[,i] <- apply(ref.dt[, celltypes == i], 1, mean)
    }
    
    df.dat = read.table(bulk.file, header =T)
    df <- df.dat[,-c(1:3)]
    df <- normalization(u, as.matrix(df))
    
    props = real_cell_proportions = read.table(prop.file, header = T)
    
    tmp.out <- data.frame(deconvolution.compare(df, ref, props),Normalization = u)
    result.out <- rbind(result.out,tmp.out)
  } 
  write.csv(result.out, out.file, row.names = F)
  
  R.squared.out <- NULL
  for(u in normalizations){
    method = result.out$method
    u.out <- NULL
    for(i in unique(method)){
      tmp.dt <- result.out[(result.out$method ==i) & (result.out$Normalization ==u),]
   
      R2 <- cor(tmp.dt$Predict.proportion, tmp.dt$True.proportion)^2
     

      RMSE <- mean((tmp.dt$Predict.proportion - tmp.dt$True.proportion)^2)
      #Rank_RMSE <- rank(RMSE,ties.method = "max")
      
      JSD <- calculate_jsd(tmp.dt$Predict.proportion,tmp.dt$True.proportion)
      
      tmp.R2 <- data.frame(R2,RMSE,JSD, method = i, Normalization = u)
      u.out <- rbind(u.out, tmp.R2) 
  
    }
   u.out$Rank_R2 <- rank(-u.out$R2,ties.method = "max")
   u.out$Rank_RMSE <- rank(u.out$RMSE,ties.method = "max")
   u.out$Rank_JSD <- rank(u.out$JSD,ties.method = "max")
   
   R.squared.out <- rbind(R.squared.out,u.out)
    
  }
 
  write.csv(R.squared.out, paste(out.file,".R2.csv",sep = ""), row.names = F)
  
}


ref.file <- "ref.file"
bulk.file <- "bulk.file"
prop.file <- "prop.file"
out.file <- "estimated.proporition.out.csv"


#normalizations <- c('none', 'z_score', 'min_max', 'col_z_score', 'col_min_max', 'QN', 'lognorm')
normalizations <- c('none','lognorm')

region.proportion.estimation(ref.file, bulk.file, prop.file, out.file, normalizations)

```

# Demo: estimating the immune cell fractions from real RRBS data

Since real RRBS data may not contain all the regions in ref.file, we recommend users to delete the regions not occurs in real RRBS data before estimating the immune cell fractions.
We apply the immune cell fraction estimation method in EpiDISH-RPC as an example.

library(EpiDISH)

    ref.file <- "ref.file"
    bulk.file <- "bulk.file"
    out.file <- "est.prop"
    
    ref.dat = read.table(ref.file, header =T)
    ref.dt <- ref.dat[,-c(1:3)]
  
    celltypes <- colnames(ref.dt) <- read.table(text =colnames(ref.dt), sep = ".")[,2]
      
    ref <- matrix(NA, ncol = length(unique(celltypes)), nrow = nrow(ref.dt))
    colnames(ref) <- unique(celltypes)
    for(i in colnames(ref)){
        ref[,i] <- apply(ref.dt[, celltypes == i], 1, mean)
    }
      
    df.dat = read.table(bulk.file, header =T, check.names = F)
    df <- df.dat[,-c(1:3)]
  
    rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)))
    predicted_cell_proportions <- epidish(as.matrix(df), as.matrix(ref),method= 'RPC')$estF
  
  
    out <- data.frame(sample.ids =colnames(df), predicted_cell_proportions)
    write.table(out, out.file, row.names = F, sep = "\t", quote = F)



# Demo: identifying cell type specific DNA methylation regions associated with phenotypes


# Performance

# Logs


- 2024-11-03: (TODO) **Manual Update**
- 2024-10-29: **Project initialization**. The first vesiion was finished in July.
