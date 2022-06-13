# cv result aggregated ----------------------------------------------------
library(tidyverse)

cv_result <- read.csv("output/nceas_blockcv_total.csv")
head(cv_result)

cv_result <- filter(cv_result, !model %in% c("GAM-unweighted", "GLM-unweighted"))

# changing the names
cv_result$model <- ifelse(cv_result$model == "RF", "RF down-sample", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "Ranger", "RF-shallow", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "GLM", "GLM-step", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "Lasso", "GLM-lasso", cv_result$model)

cv_result$model <- ifelse(cv_result$model == "MaxEnt", "MaxEnt (default)", cv_result$model)

# number of species
length(unique(cv_result$species))

mean_cv2 <-  cv_result %>%
  group_by(model, cv) %>% 
  summarise(
    ROC_mean = mean(ROC), ROC_se = 1 * (sd(ROC) / sqrt(n())),
    PRG_mean = mean(PRG), PRG_se = 1 * (sd(PRG) / sqrt(n())),
    COR_mean = mean(COR, na.rm = TRUE), COR_se = 1 * (sd(COR, na.rm = TRUE) / sqrt(n()))
  )
mean_cv2

summary(mean_cv2)

cols2 <- c(
  "SVM" = "#A588B2",
  "Ensemble" = "#35274A",
  
  "MARS" = "#E8C520",
  "GLM-step" = "#E3AB06",
  "GLM-lasso" = "#EA5C00",
  "GAM" = "#EE3B00",
  "MaxEnt (default)" = "#BE2207",
  
  "RF-shallow" = "#28868C",
  "BRT" = "#046C9A",
  "RF down-sample" = "#0B775E"
)


mean_cv3 <- filter(mean_cv2, !model %in% c("MaxEnt-tuned", "MaxEnt-spatial-tuned", 
                                           "MaxEnt-noclamp", "MaxEnt-H2",
                                           "MaxEnt-LQ", "MaxEnt-LQP"))

ggplot(data = mean_cv3, aes(x = ROC_mean, y = COR_mean, color = model)) +
  scale_color_manual(values = cols2) +
  geom_segment(aes(x = ROC_mean, y = COR_mean - COR_se, xend = ROC_mean, yend = COR_mean + COR_se), 
               colour = "gray75", alpha = 0.8, data = mean_cv3) +
  geom_segment(aes(x = ROC_mean - ROC_se, y = COR_mean, xend = ROC_mean + ROC_se, yend = COR_mean), 
               colour = "gray75", alpha = 0.8, data = mean_cv3) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(x = ROC_mean, 
                               y = COR_mean,
                               colour = model,
                               label = model),
                           force = 5,
                           data = mean_cv3) +
  facet_wrap(~ paste(str_to_title(cv), "partitioning")) +
  labs(x = expression("AUC"["ROC"]), y = "COR") +
  theme_bw(base_line_size = 0.2) +
  theme(text = element_text(size = 12, family = "Helvetica"),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(from = 0.67, to = 0.75, by = 0.02))

ggsave("fig/models_all2.jpeg", 
       width = 8, 
       height = 4,
       dpi = 600)

## Maxent highlight
cols3 <- c(
  "SVM" = "gray90",
  "Ensemble" = "gray90", # "#CFCFCF",
  
  "MARS" = "gray90",
  "GLM-step" = "gray90",
  "GLM-lasso" = "gray90",
  "GAM" = "gray90",
  
  "MaxEnt (default)" = "#BE2207",
  "MaxEnt-noclamp" = "#35274A",
  "MaxEnt-tuned" = "#E8C520",
  "MaxEnt-spatial-tuned" = "#EA5C00",
  "MaxEnt-LQ" = "#28868C",
  
  "RF-shallow" = "gray90",
  "BRT" = "gray90",
  "RF down-sample" = "gray90"
)

ggplot(data = mean_cv2, aes(x = ROC_mean, y = COR_mean, color = model)) +
  geom_segment(aes(x = ROC_mean, 
                   y = COR_mean - COR_se,
                   xend = ROC_mean, 
                   yend = COR_mean + COR_se), 
               colour = "gray85", alpha = 0.8, data = mean_cv2) +
  geom_segment(aes(x = ROC_mean - ROC_se, 
                   y = COR_mean, 
                   xend = ROC_mean + ROC_se, 
                   yend = COR_mean), 
               colour = "gray85", alpha = 0.8, data = mean_cv2) +
  geom_point(size = 2) +
  scale_color_manual(values = cols3) +
  ggrepel::geom_text_repel(aes(x = ROC_mean, 
                               y = COR_mean, 
                               colour = model, 
                               label = model),
                           force = 15, # more important
                           max.overlaps = 15,
                           data = mean_cv2) +
  facet_wrap(~ paste(str_to_title(cv), "partitioning")) +
  labs(x = expression("AUC"["ROC"]), y = "COR") +
  theme_bw(base_line_size = 0.2) +
  theme(text = element_text(size = 12, family = "Helvetica"),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(from = 0.67, to = 0.75, by = 0.02))

ggsave("fig/models_maxent2.jpeg", 
       width = 8, 
       height = 4,
       dpi = 600)


#
# Models rank -------------------------------------------------------------
library(scmamp)

cv_result2 <- filter(cv_result, !model %in% c("MaxEnt-tuned", 
                                              "MaxEnt-spatial-tuned",
                                              "MaxEnt-noclamp",
                                              "MaxEnt-LQ"))
head(cv_result2)

auc_ranks1 <- cv_result2 %>% 
  filter(cv == "random") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = ROC) %>% 
  dplyr::select(- species) %>% 
  as.matrix() %>% 
  rankMatrix() %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  setNames("AUC") %>% 
  add_rownames(var = "model") %>% 
  mutate(cv = "Random partitioning")
auc_ranks2 <- cv_result2 %>% 
  filter(cv == "spatial") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = ROC) %>% 
  dplyr::select(- species) %>% 
  as.matrix() %>% 
  rankMatrix() %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  setNames("AUC") %>% 
  add_rownames(var = "model") %>% 
  mutate(cv = "Spatial partitioning")

auc_ranks <- bind_rows(auc_ranks1, auc_ranks2)

cor_ranks1 <- cv_result2 %>% 
  filter(cv == "random") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = COR) %>% 
  dplyr::select(- species) %>% 
  as.matrix() %>% 
  rankMatrix() %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  setNames("COR") %>% 
  add_rownames(var = "model") %>% 
  mutate(cv = "Random partitioning")
cor_ranks2 <- cv_result2 %>% 
  filter(cv == "spatial") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = COR) %>% 
  dplyr::select(- species) %>% 
  as.matrix() %>% 
  rankMatrix() %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  setNames("COR") %>% 
  add_rownames(var = "model") %>% 
  mutate(cv = "Spatial partitioning")

cor_ranks <- bind_rows(cor_ranks1, cor_ranks2)

models_ranks1 <- left_join(auc_ranks, cor_ranks, by = c("model", "cv"))

p1 <- ggplot(data = models_ranks1, aes(x = AUC, y = COR, col = model)) +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(aes(label = model), alpha = 0.6) +
  scale_color_manual(values = cols2) +
  facet_wrap(~ cv) +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_bw(base_line_size = 0.1) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        text = element_text(size = 12, family = "Helvetica"),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0.1, 0, 0.2), "cm")) +
  labs(x = expression("AUC"["ROC"]))
p1


prg_ranks1 <- cv_result2 %>% 
  filter(cv == "random") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = PRG) %>% 
  dplyr::select(- species) %>% 
  as.matrix() %>% 
  rankMatrix() %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  setNames("PRG") %>% 
  add_rownames(var = "model") %>% 
  mutate(cv = "Random partitioning")
prg_ranks2 <- cv_result2 %>% 
  filter(cv == "spatial") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = PRG) %>% 
  dplyr::select(- species) %>% 
  as.matrix() %>% 
  rankMatrix() %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  setNames("PRG") %>% 
  add_rownames(var = "model") %>% 
  mutate(cv = "Spatial partitioning")

prg_ranks <- bind_rows(prg_ranks1, prg_ranks2)

models_ranks2 <- left_join(auc_ranks, prg_ranks, by = c("model", "cv"))

p2 <- ggplot(data = models_ranks2, aes(x = AUC, y = PRG, col = model)) +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(aes(label = model), alpha = 0.6) +
  scale_color_manual(values = cols2) +
  facet_wrap(~ cv) +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_bw(base_line_size = 0.1) +
  theme(text = element_text(size = 12, family = "Helvetica"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0.05, 0.1, 0.1, 0.1), "cm")) +
  labs(x = expression("AUC"["ROC"]), 
       y = expression("AUC"["PRG"]))
p2

cowplot::plot_grid(p1, p2, ncol = 1)


ggsave("fig/models_rank.jpeg", 
       width = 8, 
       height = 7.2,
       dpi = 600)


#
# The Friedman’s Aligned Rank Test ----------------------------------------
library(scmamp)

cv_result2 <- filter(cv_result, !model %in% c("MaxEnt-tuned", "MaxEnt-spatial-tuned", 
                                              "MaxEnt-noclamp", "MaxEnt-H2",
                                              "MaxEnt-LQ", "MaxEnt-LQP"))
head(cv_result2)

################################
auc_ranks1 <- cv_result2 %>% 
  filter(cv == "random") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = ROC) %>% 
  dplyr::select(- species) %>% 
  as.matrix()

friedmanAlignedRanksTest(auc_ranks1)

post_results1 <- postHocTest(
  data = auc_ranks1,
  test = "aligned ranks",
  correct = "shaffer",
  control = NULL,
  use.rank = TRUE
)

colnames(post_results1$summary) <- gsub("F.", "F ", colnames(post_results1$summary))
colnames(post_results1$summary) <- gsub("down.sample", "down-sample", colnames(post_results1$summary))
colnames(post_results1$corrected.pval) <- gsub("F.", "F ", colnames(post_results1$summary))
rownames(post_results1$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results1$summary))
colnames(post_results1$corrected.pval) <- gsub("F.", "F ", colnames(post_results1$summary))
rownames(post_results1$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results1$summary))

tprobw1 <- post_results1$corrected.pval %>% 
  as.data.frame() %>% 
  mutate(Methods = rownames(.)) %>% 
  pivot_longer(cols = seq_len(nrow(.))) %>% 
  mutate(value = round(value, 3), 
         resampling = "Random",
         Metric = "ROC")

################################
auc_ranks2 <- cv_result2 %>% 
  filter(cv == "spatial") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = ROC) %>% 
  dplyr::select(- species) %>% 
  as.matrix()

friedmanAlignedRanksTest(auc_ranks2)

post_results2 <- postHocTest(
  data = auc_ranks2,
  test = "aligned ranks",
  correct = "shaffer",
  control = NULL,
  use.rank = TRUE
)

colnames(post_results2$summary) <- gsub("F.", "F ", colnames(post_results1$summary))
colnames(post_results2$summary) <- gsub("down.sample", "down-sample", colnames(post_results1$summary))
colnames(post_results2$corrected.pval) <- gsub("F.", "F ", colnames(post_results1$summary))
rownames(post_results2$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results1$summary))
colnames(post_results2$corrected.pval) <- gsub("F.", "F ", colnames(post_results1$summary))
rownames(post_results2$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results1$summary))

tprobw2 <- post_results2$corrected.pval %>% 
  as.data.frame() %>% 
  mutate(Methods = rownames(.)) %>% 
  pivot_longer(cols = seq_len(nrow(.))) %>% 
  mutate(value = round(value, 3), 
         resampling = "Spatial",
         Metric = "ROC")

################################
auc_ranks3 <- cv_result2 %>% 
  filter(cv == "random") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = COR) %>% 
  dplyr::select(- species) %>% 
  as.matrix()

friedmanAlignedRanksTest(auc_ranks3)

post_results3 <- postHocTest(
  data = auc_ranks3,
  test = "aligned ranks",
  correct = "shaffer",
  control = NULL,
  use.rank = TRUE
)

colnames(post_results3$summary) <- gsub("F.", "F ", colnames(post_results3$summary))
colnames(post_results3$summary) <- gsub("down.sample", "down-sample", colnames(post_results3$summary))
colnames(post_results3$corrected.pval) <- gsub("F.", "F ", colnames(post_results3$summary))
rownames(post_results3$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results3$summary))
colnames(post_results3$corrected.pval) <- gsub("F.", "F ", colnames(post_results3$summary))
rownames(post_results3$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results3$summary))

tprobw3 <- post_results3$corrected.pval %>% 
  as.data.frame() %>% 
  mutate(Methods = rownames(.)) %>% 
  pivot_longer(cols = seq_len(nrow(.))) %>% 
  mutate(value = round(value, 3), 
         resampling = "Random",
         Metric = "COR")

################################
auc_ranks4 <- cv_result2 %>% 
  filter(cv == "spatial") %>% 
  pivot_wider(id_cols = species, names_from = model, values_from = COR) %>% 
  dplyr::select(- species) %>% 
  as.matrix()

friedmanAlignedRanksTest(auc_ranks4)

post_results4 <- postHocTest(
  data = auc_ranks4,
  test = "aligned ranks",
  correct = "shaffer",
  control = NULL,
  use.rank = TRUE
)

colnames(post_results4$summary) <- gsub("F.", "F ", colnames(post_results4$summary))
colnames(post_results4$summary) <- gsub("down.sample", "down-sample", colnames(post_results4$summary))
colnames(post_results4$corrected.pval) <- gsub("F.", "F ", colnames(post_results4$summary))
rownames(post_results4$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results4$summary))
colnames(post_results4$corrected.pval) <- gsub("F.", "F ", colnames(post_results4$summary))
rownames(post_results4$corrected.pval) <- gsub("down.sample", "down-sample", colnames(post_results4$summary))

tprobw4 <- post_results4$corrected.pval %>% 
  as.data.frame() %>% 
  mutate(Methods = rownames(.)) %>% 
  pivot_longer(cols = seq_len(nrow(.))) %>% 
  mutate(value = round(value, 3), 
         resampling = "Spatial",
         Metric = "COR")

################################
# c(bottom, left, top, right)
par(mfrow=c(4,1), mai = c(0, 0.1, 0, 0))
plotRanking(post_results1$corrected.pval, 
            post_results1$summary,
            alpha = 0.05, 
            cex = 1,
            decreasing = FALSE)

plotRanking(post_results2$corrected.pval, 
            post_results2$summary,
            alpha = 0.05, 
            cex = 1,
            decreasing = FALSE)

plotRanking(post_results3$corrected.pval, 
            post_results3$summary,
            alpha = 0.05, 
            cex = 1,
            decreasing = FALSE)

plotRanking(post_results4$corrected.pval, 
            post_results4$summary,
            alpha = 0.05, 
            cex = 1,
            decreasing = FALSE)

################################

alleval <- bind_rows(tprobw1, tprobw2)

ord <- rev(colnames(post_results1$summary)[order(post_results1$summary)])
mycol <- viridis::viridis(n = 30, option = "A", direction = -1)[-c(28:30)]

ggplot(data = alleval, aes(x = Methods, y = name, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), colour = "gray80", size = 2.5) +
  labs(x = NULL, y = NULL, fill = "") +
  theme_minimal() +
  facet_wrap(~resampling) +
  # viridis::scale_fill_viridis(option = "A", direction = -1, na.value  = NA,
  #                             limits = c(0, 1)) +
  scale_fill_gradientn(colours = mycol) +
  theme(axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.title.y = element_text(margin = margin(r = 10)),
        text = element_text(size = 11, family = "Helvetica"),
        axis.ticks = element_line(),
        legend.position = "top", 
        legend.direction = "horizontal",
        legend.justification =  "center",
        legend.margin = margin(0,0,0,0)
  ) +
  # guides(fill = guide_colorbar(barwidth = 15, 
  #                              barheight = 1, 
  #                              title.position = "right", 
  #                              title.vjust = 0.1,
  #                              label.position =  "top",
  #                              title = expression("AUC"["PRG"]))) +
  guides(fill = FALSE) +
  scale_x_discrete(limits = ord) +
  scale_y_discrete(limits = ord)

#
##################################################

plotRanking(post_results1$corrected.pval, 
            post_results1$summary,
            alpha = 0.05, 
            cex = 1,
            decreasing = FALSE)

# # plot models significance
drawAlgorithmGraph(post_results1$corrected.pval,
                   mean.value = post_results1$summary,
                   alpha = 0.05,
                   font.size = 10)

# plotCD(auc_ranks1, alpha=0.05, cex=1.25)

#
# perfromance rank - top 1-3 ----------------------------------------------
library(tidyverse)

cv_result <- read.csv("output/nceas_blockcv_total.csv")
head(cv_result)

cv_result <- filter(cv_result, !model %in% c("GAM-unweighted", "GLM-unweighted"))

# changing the names
cv_result$model <- ifelse(cv_result$model == "RF", "RF down-sample", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "Ranger", "RF-shallow", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "GLM", "GLM-step", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "Lasso", "GLM-lasso", cv_result$model)
cv_result$model <- ifelse(cv_result$model == "MaxEnt", "MaxEnt (default)", cv_result$model)

cv_result <- filter(cv_result, !model %in% c("MaxEnt-tuned", "MaxEnt-spatial-tuned", 
                                             "MaxEnt-noclamp", "MaxEnt-H2",
                                             "MaxEnt-LQ", "MaxEnt-LQP"))


# ensmodels <- c(
#   "GAM",
#   "GLM-lasso",
#   "BRT",
#   "RF down-sample",
#   "MaxEnt (default)",
#   "Ensemble"
# )
# cv_result <- cv_result %>%
#   filter(model %in% ensmodels)
# head(cv_result)

lastn <- function(x, n = 3, by = "COR"){
  x1 <- unique(x[, by, drop = TRUE])
  idx <- tail(order(x1), n)
  xx <- x1[idx]
  idx2 <- which(x[, by, drop = TRUE] %in% xx)
  return(x[idx2, ])
}

topModFun <- function(evalmetric, nm){
  cv_result %>% 
    filter(cv == "spatial") %>% # *** use spatial or random ****
    group_by(species) %>% 
    nest() %>% 
    mutate(df = map(data, ~lastn(x = ., n = nm, by = evalmetric)),
           models = map(df, pluck("model"))) %>% 
    select(models) %>% 
    map(unlist) %>% 
    pluck("models") %>% 
    table() %>% 
    as.data.frame() %>% 
    setNames(c("Model", "freq")) %>% 
    mutate(prop = (freq / 171) * 100,
           percent = round(prop, 1),
           metric = evalmetric,
           position = paste("Top", nm))
}

topModels <- map2(rep(c("ROC", "PRG", "COR"), each = 3), 
                  rep(c(1,2,3), 3), 
                  topModFun) %>% 
  bind_rows()

# create labes for the facet's labeller
topModels$metlabel <- factor(topModels$metric, labels = c(
  'COR',
  '"AUC"["PRG"]',
  '"AUC"["ROC"]'
))
topModels$metlabel <- fct_relevel(topModels$metlabel, c(c('"AUC"["ROC"]',
                                                          '"AUC"["PRG"]',
                                                          'COR')))
topModels

# for(i in arrange(models_info_all2, auc_mean)$Models) print(i)
morder <- c("SVM",
            "MARS",
            "RF-shallow",
            "GAM",
            "GLM-step",
            "GLM-lasso",
            "BRT",
            "RF down-sample",
            "MaxEnt (default)",
            "Ensemble")

# plot all the metrics
ggplot(data = topModels, aes(x = Model, y = forcats::fct_rev(as.factor(position)), fill = percent)) + 
  geom_tile(color = "gray") +
  facet_wrap(vars(metlabel), nrow = 3, strip.position = "right", labeller = label_parsed) +
  geom_text(aes(label = percent, colour = percent), size = 3.5) +
  viridis::scale_fill_viridis(option = "A", direction = -1) +
  viridis::scale_colour_viridis(option = "E", direction = 1, begin = 0.2, end = 0.8) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(margin = margin(r = 10)),
        text = element_text(size = 11, family = "Helvetica")) +
  guides(fill = "none", colour = "none") +
  scale_x_discrete(limits = morder)

# ggsave(sprintf("figs/%s.jpg", "fig9"),
#        units = "in",
#        dpi = 600,
#        width = 6,
#        height = 4)


#
# maxent tuned param ------------------------------------------------------
library(wesanderson)

specieslist <- read.csv("species_list.csv")

params <- read.csv("output/maxexnt_param_total.csv") %>%
  filter(speice %in% unique(specieslist$Species)) %>% 
  mutate(fold = paste0("fold", fold))
head(params)

maxparam <- function(np){
  if (np < 10) 
    classes <- "L"
  else if (np < 15) 
    classes <- "LQ"
  else if (np < 80) 
    classes <- "LQH"
  else classes <- "LQHP"
  return(classes)
}

params_def <- read_csv("output/nceas_blockcv_nocc.csv") %>% 
  mutate(model = "MaxEnt (default)",
         reg = "betamultiplier=1",
         args = sapply(.$mean_occ, maxparam)) %>% 
  dplyr::select(-mean_occ)
head(params_def)
nrow(params_def)

paramsall <- bind_rows(params, params_def) %>% 
  filter(model != "MaxEnt-tuned")
table(paramsall$model)

ggplot(data = paramsall, aes(x = args, fill = reg)) +
  geom_bar() +
  # theme_bw() +
  facet_wrap(~ model) +
  labs(fill = "Regularisation", x = "MaxEnt features", y = "Relative number") +
  scale_fill_manual(values = rev(wes_palette("Royal2"))) +
  # scale_fill_manual(values = rev(ghibli::ghibli_palette("PonyoMedium"))) +
  theme(legend.position = "top",
        text = element_text(size = 12, family = "Helvetica")) +
  scale_y_continuous(name = "Percentage",
                     labels = function(y) paste0(y * 100 / (171 * 5), "%"),
                     breaks = c(0.1, 0.2, 0.3, 0.4, 0.5) * 171 * 5) +
  labs(fill = "")

ggsave("fig/maxent_param.jpeg", 
       width = 8, 
       height = 3.5,
       dpi = 600)

# The end -----------------------------------------------------------------
