#####################################################
####### Data preparation - Cultivar Stability ####### 
#####################################################

####### Author: Vinicius Garnica
####### Date: Aug, 2024

### Load packages ------------------------------------------------------------------------------------ 

pacman::p_load(lubridate,  # must have pacman R package installed
               agricolae,
               minpack.lm,
               data.table,
               moments,
               metrica,
               skimr,
               ggpubr,
               MASS,
               tidyverse)


### Set the path for workflow

rm(list = ls())

 ### Load data sets
load("data/raw_stability.RData")

### Glimpse on data
glimpse(data)


### Functions ------------------------------------------------------------------------------------ 
### Two parameter growth models
### Monomolecular
M.2 = function(x,y0,r) {
  1 - (1 - y0)*exp(-r * x)
}
### Logistic
L.2 = function(x,y0,r) {
  (1 + ((1 - y0)/y0)*exp(-r * x))^-1
}
### Gompertz
#exp(log(y0)*exp(-r * x))
G.2 = function(x,y0,r) {
  y0 * exp( log(1/y0) * (1 - exp(-r * x)))
}

### Model fitting
fit.non.linear.2 = function(time,intensity,data=NULL,starting_par = list(y0 = 0.01, r = 0.04),maxiter = 100){ # may need to change the initial values otherwise models won't fit
  time =eval(substitute(time),data, parent.frame())
  intensity =eval(substitute(intensity),data, parent.frame())
  dt.epi = data.frame(time=time, intensity=intensity)
  
  fit_M.2 = nlsLM(intensity ~ M.2(x = time,y0=y0,r = r), 
                  start = starting_par,
                  data = dt.epi,
                  lower = c(-Inf,0), 
                  upper = c(1, Inf), 
                  control = nls.lm.control(maxiter = maxiter))
  fit_L.2 = nlsLM(intensity ~ L.2(x = time,y0=y0,r = r), 
                  start = starting_par,
                  data = dt.epi,
                  lower = c(-Inf,0), 
                  upper = c(1, Inf), 
                  control = nls.lm.control(maxiter = maxiter))
  fit_G.2 = nlsLM(intensity ~ G.2(x = time,y0=y0,r = r), 
                  start = starting_par,
                  data = dt.epi,
                  lower = c(-Inf,0), 
                  upper = c(1, Inf), 
                  control = nls.lm.control(maxiter = maxiter))
  models = data.frame(model = c("Monomolecular","Logistic","Gompertz"))  %>% 
    dplyr::mutate(y0 = dplyr::case_when(model =="Monomolecular" ~ (summary(fit_M.2)$parameters[1,1]),
                                        model =="Logistic" ~ (summary(fit_L.2)$parameters[1,1]),
                                        model =="Gompertz" ~ (summary(fit_G.2)$parameters[1,1])),
                  r =  dplyr::case_when(model =="Monomolecular" ~ (summary(fit_M.2)$parameters[2,1]),
                                        model =="Logistic" ~ (summary(fit_L.2)$parameters[2,1]),
                                        model =="Gompertz" ~ (summary(fit_G.2)$parameters[2,1])))
  fitted = bind_rows(list(
    data.frame(
      model = "Monomolecular",
      intensity = dt.epi$intensity,
      time = dt.epi$time,
      pred = predict(fit_M.2,dt.epi$time)),
    data.frame(
      model = "Logistic",
      intensity = dt.epi$intensity,
      time = dt.epi$time,
      pred = predict(fit_L.2,dt.epi$time)),
    data.frame(
      model = "Gompertz",
      intensity = dt.epi$intensity,
      time = dt.epi$time,
      pred = predict(fit_G.2,dt.epi$time))))  
  
  performance = fitted %>%
    dplyr::group_by(model) %>% tidyr::nest() %>% dplyr::ungroup() %>%
    mutate(performance = 
             map(data, ~metrica::metrics_summary(data=., 
                                                 obs = intensity, 
                                                 pred = pred,
                                                 type = "regression")),
           metrics = map(performance, ~ filter(.,Metric=="RMSE"))) %>%
    unnest(metrics) %>%
    select(-data,-performance) 
  parameters = suppressMessages(tibble(left_join(models,performance)) %>% dplyr::arrange(Score) %>% dplyr::mutate(best_model = 1:3))
  return(list(Parameters =parameters,Data=fitted))
}


### Disease measures ------------------------------------------------------------------------------------ 
dis.data = data %>% 
  select(SITE,BLOCK,CULT,PLOT,TRT,I_DATE,matches("inc"),matches("sev"),matches("date"),matches("feekes")) %>%
  mutate(I_DOY = lubridate::yday(I_DATE),
         DOY_0 = lubridate::yday(date_0),
         DOY_1 = lubridate::yday(date_1),
         DOY_2 = lubridate::yday(date_2),
         DOY_3 = lubridate::yday(date_3),
         DOY_4 = lubridate::yday(date_4),
         DOY_5 = lubridate::yday(date_5),
         DOY_6 = lubridate::yday(date_6),
         EU = factor(SITE):factor(PLOT)) %>%
  select(-matches("date")) %>% ungroup



### 1) Final percent diseased leaf area (sev) ---------------------------------------------------------------
FPDLA = dis.data %>%
  select(-matches("inc"),-matches("DOY")) %>%
  pivot_longer(cols = -c(SITE,BLOCK,CULT,PLOT,TRT,EU),names_to = c(".value","TIME"),  names_sep="_" ) %>%
  drop_na("sev") %>%
  group_by(SITE,BLOCK,CULT,PLOT,TRT,EU) %>%
  filter(TIME == max(TIME)) %>%
  ungroup() %>% # get the highest inc for each plot
  select(EU,sev)


### 2) Normalized area under the disease progress stairs (rAUDPS) (Fry, 1978) ---------------------------------------
AUDPS = dis.data %>%
  select(-matches("inc")) %>%
  pivot_longer(cols = -c(SITE,BLOCK,CULT,PLOT,TRT,I_DOY,EU),names_to = c(".value","TIME"),  names_sep="_" ) %>%
  drop_na("sev") %>%
  arrange(EU,TIME) %>% 
  group_by(SITE,BLOCK,CULT,PLOT,TRT,EU) %>%
  dplyr::summarise(audps = audps(sev,DOY,type = "absolute"),
                   first_DOY = min(DOY),
                   last_DOY = max(DOY),
                   I_DOY = max(I_DOY)) %>%
  ungroup() %>%
  mutate(raudps = audps/(last_DOY-first_DOY))

nAUDPS = AUDPS %>%
  rowwise() %>%
  mutate(naudps = round(raudps/(last_DOY-I_DOY),3)) %>%
  select(EU,audps,raudps,naudps)


### 3) Weighted Mean Absolute Rate of Disease Increase for Incidence (WMAR_INC) ---------------------------------------------
inc.long = dis.data %>%
  select(-matches("sev"),-I_DOY) %>%
  pivot_longer(cols = -c(SITE,BLOCK,CULT,PLOT,TRT,EU),names_to = c(".value","TIME"),  names_sep="_" ) %>%
  drop_na("inc") %>%
  mutate(inc = (inc/5)) %>% # create individual experimental units
  arrange(SITE,BLOCK,CULT,TIME) %>%
  select(-feekes)


disease.fit = inc.long %>%
  ungroup() %>%
  group_by(EU) %>%
  nest() %>%
  mutate(fit = 
           map(data,~fit.non.linear.2(time = DOY,intensity = inc,data=.x)))


WMAR_INC  = unnest_wider(disease.fit,fit) %>% 
  unnest(Parameters) %>%
  filter(best_model == 1) %>%
  mutate(omega = case_when(model == "Monomolecular"~ r/(2*0 + 2),
                         model == "Logistic" ~r/(2*2 + 2),
                         model == "Gompertz" ~r/(2*1 + 2)),
         omega =round(omega,3)) %>%
  select(-data, -Data,-Score,-Metric,-best_model) 

### Visually Evaluate model fit
#best=unnest_wider(disease.fit,fit) %>% 
#  unnest(Parameters) %>%
#  filter(best_model == 1) %>%
#  select(EU,model)

#left_join(best,unnest_wider(disease.fit,fit) %>% 
#            unnest(Data)) %>% filter(grepl("^PY", EU))%>%
#  ggplot(aes(x=time,y=intensity))+
#  geom_line(aes(y=pred))+
#  geom_point()+
#  facet_wrap(~EU+model)

### 4) Time Required for Incidence to reach 50% (TI_50)  ---------------------------------------------------------------

T_50 = WMAR_INC %>% 
  group_by(EU) %>%
  mutate(t50 = case_when(model == "Monomolecular" ~ -(1/r) * log((1-0.5)/(1-y0)), # boa
                         model == "Logistic" ~  -(1/r) * log(((1/0.5)-1)*(y0/(1-y0))), # boa
                         model == "Gompertz" ~ -(1/r) * log(1 - log(0.5/y0) / log(1/y0)),
                         ),
         t50 = round(t50,0)) %>%
  select(-r,-omega,-y0)


### Export data -----------------------------------------------------------------------------------------------------------

FPDLA # Final percent diseased leaf area
nAUDPS # Normalized area under the disease progress stairs
WMAR_INC # Weighted mean absolute rate of disease increase for incidence
T_50 # Time required for Incidence to reach 50%


stability_final = data %>%
  mutate(EU = factor(SITE):factor(PLOT)) %>%
  select(SITE,BLOCK,CULT,PLOT,TRT,EU) %>%
  left_join(FPDLA)%>%
  left_join(nAUDPS) %>%
  left_join(WMAR_INC) %>%
  left_join(T_50) %>%
  select(-model,-y0,-r) %>%
  mutate(omega = omega*100) # transform omega to match the 

stability_final = stability_final %>% 
  mutate_at(vars(SITE,BLOCK,CULT,PLOT,TRT,EU), as.factor) %>%
  mutate_at(vars(sev,audps,raudps,naudps,omega,t50), as.double) %>%
  mutate_at(vars(audps,raudps,omega),round,2) %>%
  data.frame() 


### Skim ------------------------------------------------------------------------------------ 
skim(stability_final)

### Save  ------------------------------------------------------------------------------------- 
save(stability_final, file = "data/stability_final.RData")
