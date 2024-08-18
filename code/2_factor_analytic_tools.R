#####################################
####### Factor Analytic Tools ####### 
#####################################

####### Author: Vinicius Garnica
####### Date: Aug, 2024

####### Some of the codes here were adapted from:
# Smith,  et al (2015). Factor analytic mixed models for the provision of grower information from national crop variety testing programs. Theor Appl Genet 128, 55–72.
# Chaves et al (2023). Employing factor analytic tools for selecting high-performance and stable tropical maize hybrids. Crop Science, 63, 1114–1125.


### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(asreml,       # must have pacman and asreml R packages installed 
               data.table,
               ggrepel,
               tidyverse)


rm(list = ls())

### Load data sets
load("data/res_stability.RData")
load("data/stability_final.RData")

### Model performance ----------------------------------------------------------------------------------------------------
res_stability$res_sev$modsel$var="sev"
res_stability$res_t50$modsel$var="t50"
res_stability$res_raudps$modsel$var="raudps"
res_stability$res_omega$modsel$var="omega"

model=rbind(res_stability$res_omega$modsel,
            res_stability$res_sev$modsel,
            res_stability$res_t50$modsel,
            res_stability$res_raudps$modsel)

write.csv(model,"results/model_selection.csv",na = "")

### Functions ----------------------------------------------------------------------------------------------------

fa.out = function(model){
  data=model$mf
  num.cultivar = nlevels(data$CULT)
  num.sites = nlevels(data$SITE)
  name.cultivar = levels(data$CULT)
  name.sites = levels(data$SITE)

  sum_model = summary(model)
  load = sum_model$varcomp[grep('fa\\d+', rownames(sum_model$varcomp)),]
  psi = sum_model$varcomp[grep('var', rownames(sum_model$varcomp)), 1]
  
  ## Original matrix of loading -------------------
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[,c('fa','component')], f = load$fa), function(x) x[,'component']))
  rownames(mat.loadings) = name.sites
  
  ###  Singular value decomposition of loadings to have  unique solution  -------------------
  svd.L = svd(mat.loadings)
  Dmat = diag(svd.L$d^2, nrow = length(svd.L$d))

  if(sum(svd.L$u[,1] < 1)/nrow(svd.L$u) > 0.5){
    svd.L$u = -1 * svd.L$u
    svd.L$v = -1 * svd.L$v
  }

  mat.loadings.star = svd.L$u

  colnames(mat.loadings.star) = colnames(mat.loadings)
  rownames(mat.loadings.star) = rownames(mat.loadings)
  
  ### Scores rotation via SVD -------------------------------------
  scores = data.frame(summary(model, coef = T)$coef.random)[
    grep('Comp', rownames(data.frame(summary(model, coef = T)$coef.random))),
  ]
  
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  scor.vec = do.call(c, lapply(split(scores, f = scores$fa), function(x) x[,'solution']))
  scor.vec.star = kronecker(sqrt(Dmat) %*% t(svd.L$v), diag(num.cultivar)) %*% scor.vec
  scor.mat.star = matrix(scor.vec.star, nrow = num.cultivar, ncol = length(unique(scores$fa)),
                         dimnames = list(name.cultivar, unique(load$fa)))
  
  ### Full variance-covariance genetic matrix and genetic correlation matrix -------------------------------------
  Gmat = mat.loadings.star%*%Dmat%*%t(mat.loadings.star) + diag(psi)
  Cmat = cov2cor(Gmat)
  
  ### Percentage of explained variance -------------------------------------
  expvar.j = matrix(ncol = ncol(mat.loadings), nrow = num.sites, 
                    dimnames = list(name.sites,colnames(mat.loadings)))
  

  for (i in 1:ncol(expvar.j)) {
    expvar.j[, i] = 100 * (Dmat[i,i] * diag(mat.loadings.star[, i] %*% t(mat.loadings.star[, i]))) / diag(Gmat)
  }
  
  if (ncol(expvar.j) > 1) {
    all = 100 * diag(mat.loadings.star %*% Dmat %*% t(mat.loadings.star)) / diag(Gmat)
    expvar.j = round(t(cbind(expvar.j, all)),1) 
  } else {
    expvar.j = round(t(expvar.j),1) 
  }
  

  ### Coefficient of variation -------------------------------------
  R = data.frame(SITE=str_sub(rownames(sum_model$varcomp[endsWith(rownames(sum_model$varcomp),"!R"),"component", drop = F]),6,-3),
                 sigma_bar = sqrt(sum_model$varcomp[endsWith(rownames(sum_model$varcomp),'!R'),1]))
  
  CV = left_join(R,(data %>% group_by(SITE) %>% summarise(x_bar = mean(response ,na.rm=TRUE)))) %>%
    mutate(CV=sigma_bar/x_bar)
  
  ### eBLUPs -------------------------------------
  blups = predict(model, classify = "SITE:CULT",maxitr=1)$pvals
  blups = blups[,-5]

  temp = data.frame("SITE" = rep(levels(data[, SITE]), each = nlevels(data[, CULT])),
                    "CULT" = rep(levels(data[, CULT]), times = nlevels(data[, SITE])),
                    "marginal" = kronecker(mat.loadings.star, diag(num.cultivar)) %*% scor.vec.star) 
  
  blups = left_join(blups, temp, by = c("SITE","CULT"))
  colnames(blups)[which(colnames(blups) == 'predicted.value')] = 'conditional'
  
  ### Overall Performance and Root Mean Square Deviation ---------------------------
  blups$L_1 = rep(mat.loadings.star[,1], each = num.cultivar)
  
  OPtest = data.frame("SITE" = rep(levels(data[, SITE]), each = nlevels(data[, CULT])),
                      "CULT" = rep(levels(data[, CULT]), times = nlevels(data[, SITE])),
                      "value"= kronecker(mat.loadings.star[,1],as.matrix(scor.mat.star[,1])))
  
  OPA = OPtest %>% group_by(CULT) %>% summarise(OP = mean(value))
  rest = (blups$marginal -(kronecker(mat.loadings.star[,1],diag(num.cultivar))) %*% scor.mat.star[,1])^2
  STA = data.frame("CULT" = rep(levels(data[, CULT]), times = nlevels(data[, SITE])),
                   "ST" = rest)
  STA = STA %>% group_by(CULT) %>% summarise(ST = sqrt(mean(ST)))
  
return(list(Gmat=Gmat, Cmat=Cmat,R=R,CV=CV,expvar.j= expvar.j,loadings=mat.loadings.star,scores=scor.mat.star,blups=blups,OPA=OPA,STA=STA))

}

### Factor Analytic Outputs  ----------------------------------------------------------------------------------------------------

fa_out_raudps = fa.out(res_stability$res_raudps$model)
fa_out_sev = fa.out(res_stability$res_sev$model)
fa_out_t50 = fa.out(res_stability$res_t50$model)
fa_out_omega = fa.out(res_stability$res_omega$model)

### Save results --------------------------------------------------------------------------------

res_fa_tools = list("omega"= fa_out_omega,
                     "sev"= fa_out_sev,
                     "t50"= fa_out_t50,
                     "raudps"=fa_out_raudps) 

save(res_fa_tools, file = "data/res_fa_tools.RData")

### Loadings -------------------------------------------------------------------------------------------
raudps = fa_out_raudps$loadings 
omega = fa_out_omega$loadings
t50 = fa_out_t50$loadings
sev = fa_out_sev$loadings

loadings = list("t50"= t50,
                "raudps"= raudps,
                "omega"= omega,
                "sev" =sev) 

save(loadings, file = "results/loadings.RData")

write.csv(reshape2::melt(loadings),"results/loadings.csv",na = "")
