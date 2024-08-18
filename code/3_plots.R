################################################
####### GxE analysis - Results and Plots ####### 
################################################

####### Author: Vinicius Garnica
####### Date: Aug, 2024

### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(asreml,       # must have pacman and asreml R packages installed 
               data.table,
               patchwork,
               viridis,
               skimr,
               ggh4x,
               biscale,
               cowplot,
               tidyverse)


rm(list = ls())

### Load data sets
load("data/res_fa_tools.RData")
load("data/stability_final.RData")


### Functions -------------------------------------------------------------------------------------------
mean_function = function(correlation_matrix) {
  correlation_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "Row") %>%
    pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%
    filter(Value != 1)
}

type_b_plot = function(correlation_matrix, var_name) {
  heatmap_df = as.data.frame(correlation_matrix)
  heatmap_df$Row = rownames(heatmap_df)
  heatmap_df_long = reshape2::melt(heatmap_df, id.vars = "Row", variable.name = "Column", value.name = "Value")
  
  eclid=hclust(dist(1-correlation_matrix),method = "complete")
  ordered_labels=eclid$labels[eclid$order]
  
  heatmap_df_long=heatmap_df_long %>%
    mutate(Row = factor(Row, levels = ordered_labels),
           Column = factor(Column, levels = ordered_labels))
  
  heatmap_df_long$VAR=var_name
  
  ggplot(heatmap_df_long, aes(x = Row, y = Column, fill = Value, label = round(Value, 2))) +
    geom_tile(color = "white") +
    geom_tile(data = filter(heatmap_df_long, Row == Column), fill = "white") +
    scale_fill_gradientn(colors = c("#075AFF", "#FFFFCC", "#FF0000"),
                         values = scales::rescale(c(-1, 1)),
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         limits = c(-1, 1),
                         na.value = "white") +
    geom_text(color = "white", size = 2.5, fontface = "bold") +
    facet_grid(cols = vars(VAR), labeller = label_parsed)+
    theme_bw() +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          legend.position = "none",
          strip.background = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 15),
          strip.text = element_text(size = 15, hjust = 0)) +
    labs(fill = "Correlation")
}


optimize_join = function(data, var_name) {
  data$OP %>%
    left_join(data$ST, by = "CULT") %>%
    mutate(VAR = factor(var_name))
}

### Data summaries -------------------------------------------------------------------------------------------

data_summary = stability_final %>%
  group_by(SITE) %>%
  select(sev, raudps, omega, t50) %>%
  skim() %>%
  tibble::as_tibble()

write.csv(data_summary, "results/data_summary.csv", na = "")

### Coefficient of Variation -------------------------------------------------------------------------------------------------------------

cv_data = bind_rows(
  res_fa_tools$omega$CV %>% mutate(var = "omega"),
  res_fa_tools$raudps$CV %>% mutate(var = "raudps"),
  res_fa_tools$t50$CV %>% mutate(var = "t50"),
  res_fa_tools$sev$CV %>% mutate(var = "sev")
)

write.csv(cv_data, "results/CV.csv", na = "")

### Explained Variation -------------------------------------------------------------------------------------------------------------

expvar_data = bind_rows(
  res_fa_tools$omega$expvar.j %>% as.data.frame() %>% mutate(var = "omega"),
  res_fa_tools$raudps$expvar.j %>% as.data.frame() %>% mutate(var = "raudps"),
  res_fa_tools$t50$expvar.j %>% as.data.frame() %>% mutate(var = "t50"),
  res_fa_tools$sev$expvar.j %>% as.data.frame() %>% mutate(var = "sev")
)

write.csv(expvar_data, "results/exp_var.csv", na = "")



### Type B cultivar correlations-------------------------------------------------------------------------------------------

mean_function(res_fa_tools$omega$Cmat) %>%
  summarise(mean(Value), min(Value), max(Value))

plots_list = list(
  type_b_plot(res_fa_tools$raudps$Cmat, "rAUDPS"),
  type_b_plot(res_fa_tools$sev$Cmat, "SEV"),
  type_b_plot(res_fa_tools$t50$Cmat, "T[50]"),
  type_b_plot(res_fa_tools$omega$Cmat, "omega")
)


plot_grid(
  wrap_plots(plots_list, ncol = 2)
)

ggsave("results/plots/fig2.tiff",width = 11,height = 11)



### Latent regressions --------------------------------------------------------------------------------------
# Data wrangling
dt_plot1 = rbind(
            optimize_join(res_fa_tools$raudps, "rAUDPS"),
            optimize_join(res_fa_tools$sev, "SEV"),
            optimize_join(res_fa_tools$omega, "omega"),
            optimize_join(res_fa_tools$t50, "T[50]"))


# Top 5 Cultivars in OP and Stability by Disease Metric
dt_plot1 %>% 
  group_by(VAR) %>%
  arrange(desc(ST)) %>%
  slice_head(n = 5)

dt_plot1 %>% 
  group_by(VAR) %>%
  arrange(ST) %>%
  slice_head(n = 5)


# Prepare Data for OP versus RMSD Plot
dt_plot1 = dt_plot1 %>%
  mutate(VAR = factor(VAR, levels = c("rAUDPS", "SEV", "T[50]", "omega")))

write.csv(dt_plot1, "results/op_rsmd.csv", na = "")

p1 = ggplot(dt_plot1, aes(x = ST, y = OP, color = CULT)) +
  geom_point(size = 2) +
  geom_text(data = filter(dt_plot1, CULT %in% c("CP9606", "DG Shirley", "USG 3230", "Jamestown")), 
            aes(label = CULT), size = 3.5, vjust = -0.5) +
  scale_color_manual(values = c("CP9606" = "blue", "DG Shirley" = "black", "USG 3230" = "green", "Jamestown" = "purple")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(clip = "off") +
  xlab(expression(italic("RMSD"))) +
  ylab(expression(italic("OP"))) +
  facet_wrap(~VAR, scales = "free", ncol = 1, labeller = label_parsed) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.title.x = element_text(size = 15, margin = margin(t = 15)),
    axis.title.y = element_text(size = 15, margin = margin(l = 15)),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    panel.grid = element_line(color = "grey50", linewidth = 0.1),
    axis.ticks.length.y = unit(0, "pt"),
    strip.text = element_text(size = 15, hjust = 0, colour = "black"),
    legend.position = "none"
  )

# EBLUP versus FA1
res_fa_tools$raudps$blups$VAR = as.factor("rAUDPS")
res_fa_tools$sev$blups$VAR = as.factor("SEV")
res_fa_tools$t50$blups$VAR = as.factor("T[50]")
res_fa_tools$omega$blups$VAR = as.factor("omega")

# Prepare Data for EBLUP versus FA1
dt_plot2 = rbind(
 res_fa_tools$raudps$blups,
  res_fa_tools$sev$blups,
  res_fa_tools$t50$blups,
  res_fa_tools$omega$blups
)

dt_plot2$VAR = factor(dt_plot2$VAR,levels = c(
  "rAUDPS", 
  "SEV", 
  "T[50]",
  "omega"))

write.csv(dt_plot1, "results/eblup_L1.csv", na = "")

# Plot EBLUP versus FA.1
p2 = dt_plot2 %>%
  group_by(VAR) %>%
  mutate(mean_loading = mean(L_1)) %>%
  filter(CULT %in% c("CP9606", "DG Shirley", "USG 3230", "Jamestown")) %>%
  ggplot(aes(x = L_1, y = marginal, color = CULT, group = CULT)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, fill = NA, linewidth = 0.5) +
  scale_color_manual(values = c("CP9606" = "blue", "DG Shirley" = "black", "USG 3230" = "green", "Jamestown" = "purple")) +
  guides(color = guide_legend(title.position = "top", title = "Cultivar"),
         shape = guide_legend(title.position = "top", title = "Region")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = expression(paste("EBLUP")))) +
  facet_wrap(~factor(VAR), scales = "free", ncol = 1, labeller = label_parsed) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(aes(xintercept = mean_loading), linetype = "dashed") +
  xlab(expression(paste("Rotated ", hat('λ')[1]))) +
  ylab(NULL) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.title.x = element_text(size = 15, margin = margin(t = 8)),
    axis.title.y = element_blank(),
    axis.text.y.right = element_text(),
    axis.title.y.right = element_text(size = 15, margin = margin(l = 10)),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.37),
    panel.grid = element_line(color = "grey50", linewidth = 0.1),
    axis.ticks.length.y = unit(0, "pt"),
    strip.text = element_text(size = 15, hjust = 0, colour = "white"),
    legend.position = "none"
  )

plots_list = list(
  p_VAR1 = p1,
  p_VAR2 = p2)

plot_grid(
  wrap_plots(plots_list, ncol = 2)
)

ggsave("results/plots/fig3.tiff",width =7,height = 10,units = "in",limitsize = FALSE)


### Host resistance ------------------------------------------------------------------------------------ 
# Classify host resistance based on quartiles of final disease severity

res_fa_tools$sev$OPA= res_fa_tools$sev$OPA %>%
   mutate(Percentile = ntile(OP, 3),
          resistance=as.factor(case_when(Percentile==1~'MR',
                                         Percentile==2~'MS',
                                         Percentile==3~'S')))

### Resistance class Stability ------------------------------------------------------------------------------------ 

data=left_join(dt_plot2 %>%
                 filter(VAR=="SEV"),
               res_fa_tools$sev$OPA %>% 
                 select(-OP,-Percentile) %>% arrange(resistance)) 

### Calculate GEI
Cmat = res_fa_tools$sev$Cmat
diag(Cmat) = NA
aux = data.frame(GEI = rowMeans(Cmat, na.rm = TRUE), SITE = row.names(Cmat))

data = data %>%
  left_join(aux, by = "SITE") %>%
  group_by(SITE) %>%
  mutate(ranking = rank(marginal)) %>%
  group_by(resistance, SITE) %>%
  mutate(class_rank = mean(ranking, na.rm = TRUE))

### Disease Susceptibility Rankings Plot ------------------------------------------------------------------------------------ 

data %>%
  ggplot(aes(x = reorder(SITE, -GEI))) +
  geom_line(aes(y = class_rank, group = resistance, linetype = resistance), linewidth = 1.2) +
  geom_line(aes(y = ranking, group = CULT, linetype = resistance), alpha = 0.08) +
  theme_bw() +
  scale_y_reverse() +
  scale_linetype_manual(values = c("solid","dotdash","dotted" )) +  
  labs(
    y = "EBLUP rank",
    x = "Environment (lower to higher GEI)",
    linetype = "Cultivar reaction"
  ) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, "lines"),  
    legend.text = element_text(size = 11),  
    legend.spacing.y = unit(0.5, "cm"), 
    legend.title = element_text(size = 13)  
  )

ggsave("results/plots/fig4.tiff",width =7.,height = 4.5,units = "in",limitsize = FALSE)

### Overall rankings

dat = bi_class(data %>% 
                  group_by(SITE) %>%
                  mutate(ranking = rank(desc(marginal))), x = L_1, y =ranking  , style = "quantile", dim = 3)

main = ggplot(dat) +
  geom_tile(aes(reorder(SITE, -L_1), reorder(CULT, ranking), group = CULT, fill =bi_class)) +
  bi_scale_fill(pal = "DkViolet2", dim =3)+
  theme_bw() +
  guides(fill="none",alpha="none")+
  labs(y = "Cultivar", x = "Environment", fill = "Rank") +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


legend = bi_legend(pal = "DkViolet2",
                    dim = 3,
                    xlab = "Higher first loading",
                    ylab = "Higher ranks",
                    size = 10)

main + legend + plot_layout(ncol =  2,widths =  c(1,0.3))

ggsave("results/plots/fig5.tiff",width =7,height = 4.5,units = "in",limitsize = FALSE)


### Ranking by cultivar

a=data[data$CULT=="CP9606",]
colSums(table(a$SITE,a$ranking))

