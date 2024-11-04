##########################################################
####### GxE analysis - Model fitting and selection ####### 
##########################################################

####### Author: Vinicius Garnica
####### Date: Aug, 2024

### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(asreml,       # must have pacman and asreml R packages installed 
               data.table,
               skimr,
               tidyverse)

### Set the path for workflow
rm(list = ls())

### Load data sets
load("data/stability_final.RData")

### Rename and glimpse data
data = stability_final 
glimpse(data)
write.csv(data,"data/stability_final.csv")

### Number of cultivars and environments  ------------------------------------------------------------------------------------- 
num.cultivar = nlevels(data$CULT)
num.sites = nlevels(data$SITE)

name.cultivar = levels(data$CULT)
name.sites = levels(data$SITE)


### Model fitting function ----------------------------------------------------------------------------------------------------
model_fit = function(response,data=NULL){
    response =eval(substitute(data$response),data, parent.frame())
    BLOCK =eval(substitute(data$BLOCK),data, parent.frame())
    CULT =eval(substitute(data$CULT),data, parent.frame())
    SITE =eval(substitute(data$SITE),data, parent.frame())
    data = data.frame(data=data, response=response, CULT=CULT,SITE=SITE,BLOCK=BLOCK)
    
## Model 1: (G:CS/R:IDV)  --------------------
mm1 = asreml(fixed = response ~ SITE + SITE:BLOCK, 
             random = ~ CULT + SITE:CULT,
             data = data, na.action = na.method(x="include", y = "include"), maxit = 100)

mm1 = update.asreml(mm1)
sum1 = summary(mm1)$varcomp;sum1
aic1 = summary(mm1)$aic;aic1
logl1 = summary(mm1)$loglik;logl1

predm1_sed = predict(mm1, classify = "CULT", sed = T)$avsed


## Model 2: (G:CS/R:IDH)  --------------------
mm2 = asreml(fixed = response ~ SITE + SITE:BLOCK, 
             random = ~ CULT + SITE:CULT,
             residual=~dsum(~id(units)|SITE),
             data = data, na.action = na.method(x="include", y = "include"), maxit = 100)

mm2 = update.asreml(mm2)
sum2 = summary(mm2)$varcomp;sum2
aic2 = summary(mm2)$aic;aic2
logl2 = summary(mm2)$loglik;logl2

predm2_sed = predict(mm2, classify = "CULT", sed = T)$avsed

  
## Model 3: (G:CSH/R:IDH)    --------------------       
mm3 = asreml(fixed = response ~ SITE  + SITE:BLOCK, 
             random = ~ corh(SITE):CULT,
             residual=~dsum(~id(units)|SITE),
             data = data, na.action = na.method(x="include", y = "include"), maxit = 100)

mm3=update.asreml(mm3)
mm3=update.asreml(mm3)

sum3 = summary(mm3)$varcomp;sum3
aic3 = summary(mm3)$aic;aic3
logl3 = summary(mm3)$loglik;logl3

predm3_sed = predict(mm3, classify = "CULT")$avsed


## Model 4: (G:US/R:IDH)    --------------------  US matrix is not positive definite. Only runs on 17 environments
#mm4 = asreml(fixed = response ~ SITE+ SITE:BLOCK, 
#             random = ~ us(SITE):CULT,
#             residual=~dsum(~id(units)|SITE),
#             data = b, na.action = na.method(x="include", y = "include"), maxit = 100,ai.sing = TRUE)

#mm4=update.asreml(mm4)
#mm4=update.asreml(mm4)

#sum4 = summary(mm4)$varcomp;sum4
#aic4 = summary(mm4)$aic;aic4
#logl4 = summary(mm4)$loglik;logl4

#Gvcov = matrix(0, num.sites, num.sites)
#Gvcov[upper.tri(Gvcov, diag = T)] = sum4[grep('CULT',rownames(sum4)),1]
#Gvcov[lower.tri(Gvcov, diag = F)] = t(Gvcov)[lower.tri(t(Gvcov), diag = F)]
#gencorr = cov2cor(Gvcov)

#predm4_sed = predict(mm4, classify = "CULT")$avsed

## Model 5: (G:FA1/R:IDH)  --------------------
mm5 = asreml(fixed = response ~ SITE + SITE:BLOCK, 
             random =~ fa(SITE,1):CULT,
             residual=~dsum(~id(units)|SITE),
             data = data, 
             na.action = na.method(x="include", y = "include"), maxit = 100,
             predict = predict.asreml(classify = "CULT:SITE"))

mm5=update.asreml(mm5)
mm5=update.asreml(mm5)
mm5=update.asreml(mm5)

sum5 = summary(mm5)$varcomp;sum5
aic5 = summary(mm5)$aic;aic5
logl5 = summary(mm5)$loglik;logl5

fa1.loadings = sum5[endsWith(rownames(sum5),'fa1'),1]
mat.loadings = as.matrix(cbind(fa1.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1")
psi = as.matrix(diag(sum5[endsWith(rownames(sum5),'var'),1])) 
lamb_lamb_t = mat.loadings.star %*% t(mat.loadings.star)
Gvcov5 = lamb_lamb_t + psi

expvar5 = round((sum(diag(lamb_lamb_t))/sum(diag(Gvcov5)))*100,1)

predm5_sed = predict(mm5, classify = "CULT")$avsed

## Model 6: (G:FA2/R:IDH)   --------------------     
mm6 = asreml(fixed = response ~ SITE+ SITE:BLOCK, 
             random =~ fa(SITE,2):CULT,
             residual=~dsum(~id(units)|SITE),
             data = data, 
             na.action = na.method(x="include", y = "include"), maxit = 100,
             predict = predict.asreml(classify = "CULT:SITE"))

mm6=update.asreml(mm6)
mm6=update.asreml(mm6)
mm6=update.asreml(mm6)
mm6=update.asreml(mm6)

sum6 = summary(mm6)$varcomp;sum6
aic6 = summary(mm6)$aic;aic6
logl6 = summary(mm6)$loglik;logl6

fa1.loadings = sum6[endsWith(rownames(sum6),'fa1'),1]
fa2.loadings = sum6[endsWith(rownames(sum6),'fa2'),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2')
psi = as.matrix(diag(sum6[endsWith(rownames(sum6),'var'),1])) 
lamb_lamb_t = mat.loadings.star %*% t(mat.loadings.star)
Gvcov6 = lamb_lamb_t + psi

expvar6 = round((sum(diag(lamb_lamb_t))/sum(diag(Gvcov6)))*100,1)

predm6_sed = predict(mm6, classify = "CULT")$avsed


## Model 7: (G:FA3/R:IDH)   --------------------      
mm7 = asreml(fixed = response~ SITE + SITE:BLOCK, 
             random =~ fa(SITE,3):CULT,
             residual=~dsum(~id(units)|SITE),
             data = data, 
             na.action = na.method(x="include", y = "include"), maxit = 100,
             predict = predict.asreml(classify = "CULT:SITE"))

mm7=update.asreml(mm7)
mm7=update.asreml(mm7)
mm7=update.asreml(mm7)
mm7=update.asreml(mm7)

sum7 = summary(mm7)$varcomp;sum7
aic7 = summary(mm7)$aic;aic7
logl7 = summary(mm7)$loglik;logl7

fa1.loadings = sum7[endsWith(rownames(sum7),'fa1'),1]
fa2.loadings = sum7[endsWith(rownames(sum7),'fa2'),1]
fa3.loadings = sum7[endsWith(rownames(sum7),'fa3'),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings, fa3.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2','fa3')
psi = as.matrix(diag(sum7[endsWith(rownames(sum7),'var'),1])) 
lamb_lamb_t = mat.loadings.star %*% t(mat.loadings.star)
Gvcov7 = lamb_lamb_t + psi

expvar7 = round((sum(diag(lamb_lamb_t))/sum(diag(Gvcov7)))*100,1)

predm7_sed = predict(mm7, classify = "CULT")$avsed


# Model Selection --------------------
modsel = data.frame('model' = as.factor(paste0('mm',seq(1:7))),
                    'aic' = c(aic1,aic2,aic3,NA,aic5,aic6,aic7),
                    'loglik' = c(logl1,logl2,logl3,NA,logl5,logl6,logl7),
                    'avsed'= c(predm1_sed[2],predm2_sed[2],predm3_sed,NA,
                               predm5_sed,predm6_sed,predm7_sed), 
                    'ExpVar' = c('-','-','-','-',expvar5,expvar6,expvar7)) %>% 
  arrange(aic) %>%                                                              # model selection criterium
  dplyr::mutate(best_model = 1:7)

models = list(mm1=mm1, mm2=mm2, mm3=mm3,mm4=NA, mm5=mm5, mm6=mm6, mm7=mm7)

best_model_name = modsel$model[1] 
best_model = models[["mm7"]] # change model


return(list(modsel=modsel,model=best_model))
}


### Fitting models ----------------------------------------------------------------------------------------------------

res_raudps = model_fit(response = "raudps",data)
res_sev = model_fit(response = "sev",data)
res_t50 = model_fit(response = "t50",data)
res_omega = model_fit(response = "omega",data)

### Save --------------------------------------------------------------------------------

res_stability = list("res_raudps"=res_raudps,
                     "res_sev"= res_sev,
                     "res_t50"= res_t50,
                     "res_omega"= res_omega) 

save(res_stability, file = "data/res_stability.RData")




