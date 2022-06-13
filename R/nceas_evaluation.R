library(tidyverse)
library(precrec)
library(prg)

fli <- list.files("output/nceas_model_output/")
# subset the incomplete folds
inc1 <- strsplit(fli, "_") %>% 
  map_chr(pluck, 2) %>% 
  table()

# this is remove model with incomplete folds
# change the 14 if you use fewer number of models
num_model <- 14 * 5 # number of complete folds: models * folds

inc2 <- names(which(inc1 < num_model))
inc3 <- which(! map_chr(strsplit(fli, "_"), pluck, 2) %in% inc2)
cat("Incomplete species:\n", inc2, "\n")
cat("Number of complete species:", length(inc3) / num_model, "\n")

cv_result2 <- map(file.path("output/nceas_model_output", fli[inc3]), read.csv, stringsAsFactors = FALSE) %>% 
  do.call(rbind.data.frame, .)
head(cv_result2)

tm <- Sys.time()
mean_cv2 <- cv_result2 %>% 
  group_by(species, model, cv) %>% 
  nest() %>% 
  mutate(COR = map_dbl(data, function(x) cor(x$score, x$label, method = "pearson")),
         ROC = map_dbl(data, function(x) precrec::auc(precrec::evalmod(scores = x$score, labels = x$label))[1,4]),
         PRG = map_dbl(data, function(x) prg::calc_auprg(prg::create_prg_curve(labels = x$label, pos_scores = x$score)))) %>% 
  dplyr::select(- data)
Sys.time() - tm

write.csv(mean_cv2, "output/nceas_blockcv_total.csv", row.names = FALSE)

sessionInfo()
