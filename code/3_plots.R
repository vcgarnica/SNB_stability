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

type_b_plot_upper_triangle = function(correlation_matrix, var_name) {
  heatmap_df = as.data.frame(correlation_matrix)
  heatmap_df$Row = rownames(heatmap_df)

  heatmap_df_long = melt(heatmap_df, id.vars = "Row", variable.name = "Column", value.name = "Value")
  

  eclid = hclust(dist(1 - correlation_matrix), method = "complete")
  ordered_labels = eclid$labels[eclid$order]
  
  heatmap_df_long = heatmap_df_long %>%
    mutate(Row = factor(Row, levels = ordered_labels),
           Column = factor(Column, levels = ordered_labels))
  

  heatmap_df_long = heatmap_df_long %>%
    filter(as.numeric(Row) < as.numeric(Column))
 
  heatmap_df_long$VAR = var_name

    ggplot(heatmap_df_long, aes(x = Row, y = Column, fill = Value, label = round(Value, 2))) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("#8c8a6c","#b59a7a","#cbb89f" ,"#d3d3d3","#8fb1c2" ,"#488fb0" ,"#311e3b"),
#   scale_fill_gradientn(colors = c("#8c8a6c","#b59a7a" ,"#d3d3d3","#8fb1c2" ,"#488fb0" ),
#   scale_fill_gradientn(colors = c( "#9e3547","#ba8890" ,"#d3d3d3" ,"#4279b0" ,"#3a4e78" ),
#      scale_fill_gradientn(colors = c( "#9e3547", "#ba8890", "#d3d3d3","#4279b0","#311e3b"),
                         values = scales::rescale(c(-1, 1)),
                         breaks = c(-1, 0,  1),
                         limits = c(-1, 1),
                         na.value = "white") +
    geom_text(color = "white", size = 2.3, fontface = "bold") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          legend.position = "right",
          axis.title = element_blank(),
          plot.title = element_text(size = 15),
          strip.text = element_text(size = 15, hjust = 0)) +
    labs(fill = expression(italic(hat(r)[g])))
}


custom_pal3 <- c(
  "1-1" = "#cbb89f", # low x, medium y (light tan)
  "2-1" = "#b59a7a", # medium x, medium y (warm beige)
  "3-1" = "#8c8a6c", # high x, medium y (muted mauve)
  
  "1-2" = "#c7c7c7", # low x, low y (subtly darker gray)
  "2-2" = "#d6b8b5", # medium x, low y (soft pink-beige with a gray tint)
  "3-2" = "#7a6b84", # high x, low y (muted rose)
  
  "1-3" = "#8fb9d6", # low x, high y (light blue)
  "2-3" = "#4279b0",# medium x, high y (soft sky blue)
  "3-3" = "#311e3b"  # high x, high y (medium blue)
)
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

mean_function(res_fa_tools$sev$Cmat) %>%
  summarise(round(mean(Value),2), round(min(Value),2), round(max(Value),2))

combined_plot = 
  type_b_plot_upper_triangle(res_fa_tools$sev$Cmat, "SEV") +theme(legend.position = "none")+
  type_b_plot_upper_triangle(res_fa_tools$raudps$Cmat, "rAUDPS") +theme(legend.position = "none")+
  type_b_plot_upper_triangle(res_fa_tools$t50$Cmat, "T[50]") +theme(legend.position = "none")+
  type_b_plot_upper_triangle(res_fa_tools$omega$Cmat, "omega") +
  plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 14))

combined_plot & plot_layout(guides = "collect")

ggsave("results/plots/fig2.tiff",width = 12.5,height = 11)


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
  mutate(VAR = factor(VAR, levels = c("SEV","rAUDPS",  "T[50]", "omega")))

write.csv(dt_plot1, "results/op_rsmd.csv", na = "")


p1 = ggplot(dt_plot1, aes(x = ST, y = OP, color = CULT)) +
  geom_point(size = 3.7,alpha=0.75) +
  geom_text(data = filter(dt_plot1, CULT %in% c("TURBO", "DG Shirley", "USG 3230", "Jamestown")), 
            aes(label = CULT), size = 4, vjust = -0.5) +
 # scale_color_manual(values = c("TURBO" = "blue", "DG Shirley" = "black", "USG 3230" = "green", "Jamestown" = "purple")) +
  scale_color_manual(values = c("TURBO" = "#7a6b84", "Jamestown" = "#8c8a6c", "USG 3230" =  "#488fb0", "DG Shirley" = "#311e3b")) +
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
 res_fa_tools$sev$blups,
 res_fa_tools$raudps$blups,
  res_fa_tools$t50$blups,
  res_fa_tools$omega$blups
)

dt_plot2$VAR = factor(dt_plot2$VAR,levels = c(
  "SEV", 
  "rAUDPS", 
  "T[50]",
  "omega"))

write.csv(dt_plot1, "results/eblup_L1.csv", na = "")

# Plot EBLUP versus FA.1
p2 = dt_plot2 %>%
  group_by(VAR) %>%
  mutate(mean_loading = mean(L_1)) %>%
  filter(CULT %in% c("TURBO", "DG Shirley", "USG 3230", "Jamestown")) %>%
  ggplot(aes(x = L_1, y = marginal, color = CULT, group = CULT)) +
  geom_point(size = 3.7,alpha=0.6) +
  geom_smooth(method = lm, fill = NA, linewidth = 0.7,alpha = 0.5) +
 # scale_color_manual(values = c("TURBO" = "blue", "DG Shirley" = "black", "USG 3230" = "green", "Jamestown" = "purple")) +
  scale_color_manual(values = c("TURBO" = "#7a6b84", "Jamestown" = "#8c8a6c", "USG 3230" =  "#488fb0", "DG Shirley" = "#311e3b")) +
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
  geom_line(aes(y = class_rank, group = resistance, color=resistance,linetype=resistance),linewidth = 1.2) +
  geom_line(aes(y = ranking, group = CULT,color=resistance,linetype=resistance),linewidth = 0.5,alpha=0.1) +
  theme_bw() +
  scale_y_reverse() +
  scale_color_manual(values = c("#8c8a6c", "#ba8890","#488fb0"))+
  labs(
    y = "EBLUP rank",
    x = "Environment (lower to higher GEI)",
    color = "Cultivar Reaction", linetype = "Cultivar Reaction"
  ) +
  theme(
    axis.title.x = element_text(size = 8, margin = margin(t = 10)),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 8, margin = margin(r = 6)),
    axis.text.y = element_text(size = 6, margin = margin(r = 6)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, "lines"),  
    legend.text = element_text(size = 7),  
    legend.spacing.y = unit(0.5, "cm"), 
    legend.title = element_text(size = 8)  
  )


ggsave("results/plots/fig4.tiff",width =5,height = 3.2,units = "in",limitsize = FALSE)

### Overall rankings

dat = bi_class(data %>% 
                  group_by(SITE) %>%
                  mutate(ranking = rank(desc(marginal))), x = L_1, y =ranking  , style = "quantile", dim = 3)


custom_pal3 <- c(
  "1-1" = "#cbb89f", # low x, medium y (light tan)
  "2-1" = "#b59a7a", # medium x, medium y (warm beige)
  "3-1" = "#8c8a6c", # high x, medium y (muted mauve)
  
  "1-2" = "#e0e0e0", # low x, low y (subtly darker gray)
  "2-2" = "#d6b8b5", # medium x, low y (soft pink-beige with a gray tint)
  "3-2" = "#7a6b84", # high x, low y (muted rose)
  
  "1-3" = "#8fb9d6", # low x, high y (light blue)
  "2-3" = "#4279b0",# medium x, high y (soft sky blue)
  "3-3" = "#311e3b"  # high x, high y (medium blue)
)

bi_pal(pal = custom_pal3, dim = 3,flip_axes = TRUE, rotate_pal = TRUE)


main = ggplot(dat) +
  geom_tile(aes(reorder(SITE,L_1), reorder(CULT, ranking), group = CULT, fill =bi_class)) +
  bi_scale_fill(pal = custom_pal3, dim =3, flip_axes = TRUE, rotate_pal = TRUE)+
  theme_bw() +
  guides(fill="none",alpha="none")+
  labs(y = "Cultivar", x = "Environment", fill = "Rank") +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


legend = bi_legend(pal = custom_pal3,flip_axes = F,
                    dim = 3, 
                    xlab = "Higher \u03BB\u2081",
                    ylab = "Higher EBLUP ranks",
                    size = 10)

main + legend + plot_layout(ncol =  2,widths =  c(1,0.3))

ggsave("results/plots/fig5.tiff",width =7,height = 4.5,units = "in",limitsize = FALSE)

### Ranking by cultivar

a=data[data$CULT=="TURBO",]
colSums(table(a$SITE,a$ranking))

