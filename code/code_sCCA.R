library(compositions)
library(PMA)
library(tidyverse)


apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
  nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
  if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
  return(x)
}

produce_cca<- function(OTU, METABO, dolog = FALSE){
  X <- compositions::clr(OTU)
  colnames(X) <- paste0("X_", 1:ncol(X))
  X <- apply_nzv_filter(X)


  Y <- unname(METABO)
  colnames(Y) <- paste0("Y_", 1:ncol(Y))
  if(dolog){
    Y <- log(Y+1) %>% apply_nzv_filter()
  } else {
    Y <- apply_nzv_filter(Y)
  }

  cca_res <- list() # return empty if any problems during cca
  set.seed(123)
  perm.out = try(PMA::CCA.permute(x = X, z = Y, typex="standard",typez="standard",
                                  nperms=25, penaltyxs=seq(0.01, 1,0.1), penaltyzs=seq(0.01, 1,0.1)))
  if(!is(perm.out, "try-error")){
    cca=PMA::CCA(x = X, z = Y, typex="standard",typez="standard", penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
                 v=perm.out$v.init,K=2)

    u = cca$u
    rownames(u) <- colnames(X)
    v = cca$v
    rownames(v) <- colnames(Y)
    cca_res <- list(u = u, v = v)
  }
  return(list(cca_res = cca_res, perm.out = perm.out, cca=cca))
}


retrieve_association_cca <- function(cca_res, index_association_X, index_association_Y){


  selected_features_X <- cca_res$u %>% as.data.frame() %>% rownames_to_column("feature") %>%
    filter(V1 != 0 | V2 != 0) %>%
    pull(feature)

  selected_features_Y <- cca_res$v %>% as.data.frame() %>% rownames_to_column("feature") %>%
    filter(V1 != 0 | V2 != 0) %>%
    pull(feature)

  colnames_index_association_X <- paste0("X_", index_association_X)
  colnames_index_association_Y <- paste0("Y_", index_association_Y)


  matching_X <- data.frame(feature = rownames(cca_res$u)) %>%
    mutate(association = ifelse(feature %in% colnames_index_association_X, 1, 0)) %>%
    mutate(selected = ifelse(feature %in% selected_features_X, 1, 0)) %>%
    mutate(block = "X")

  matching_Y <- data.frame(feature = rownames(cca_res$v)) %>%
    mutate(association = ifelse(feature %in% colnames_index_association_Y, 1, 0)) %>%
    mutate(selected = ifelse(feature %in% selected_features_Y, 1, 0)) %>%
    mutate(block = "Y")

  binded <- rbind(matching_X, matching_Y) %>% mutate(truth = association, pred = selected)

  conf_tab <- caret::confusionMatrix(data = factor(binded$pred, levels = c(0,1)),
                                     reference = factor(binded$truth, levels = c(0,1)),
                                     positive = "1")

  conf_tab_vec <- conf_tab$table %>% as.data.frame() %>% mutate(Class = case_when(
    Prediction == 0 & Reference == 0 ~ "TN",
    Prediction == 0 & Reference == 1 ~ "FN",
    Prediction == 1 & Reference == 0 ~ "FP",
    Prediction == 1 & Reference == 1 ~ "TP"
  )) %>% dplyr::select(Class, Freq) %>%
    column_to_rownames("Class") %>% t %>% as_tibble() %>% c()

  TP <- conf_tab_vec$TP
  TN <- conf_tab_vec$TN
  FP <- conf_tab_vec$FP
  FN <- conf_tab_vec$FN


  metrics <- list()
  metrics[["Sparsity"]] <- sum(binded$pred)/length(binded$pred)

  metrics[["F1"]] <- conf_tab$byClass[["F1"]]
  metrics[["Specificity"]] <- conf_tab$byClass[["Specificity"]]
  metrics[["Sensitivity"]] <- conf_tab$byClass[["Sensitivity"]]

  metrics[["MCC"]] <- try((TP*TN - FP*FN)/(
    (TP + FP)*(TP+FN)*(TN+FP)*(TN+FN)
  ))

  results <- list("conf_tab" = conf_tab, "metrics" = metrics)
  return(results)
}


# apply

data_autism <- readRDS("synthetic_data_Autism_alter_with_index_associations.RDS")
data_adenomas <- readRDS("synthetic_data_Adenomas_alter_with_index_associations.RDS")
data_konzo <- readRDS("synthetic_data_from_Konzo_dataset_alter_with_index_associations_full.RDS")

names(data_autism) <- paste0("Autism_", seq_along(data_autism))
names(data_adenomas) <- paste0("Adenomas_", seq_along(data_adenomas))
names(data_konzo) <- paste0("Konzo_", seq_along(data_konzo))

all_data <- c(data_autism, data_adenomas, data_konzo)

cca_all <- parallel::mclapply(seq_along(all_data[1]),function(i){
  print(paste0("scca: ", i));
  OTU = all_data[[i]]$Simulated.Microbiotes;
  METABO = all_data[[i]]$Simulated.Metabolites;

  produce_cca(OTU = OTU,
              METABO = METABO)
}, mc.cores = 50)

cca_conf_tab_all <- parallel::mclapply(seq_along(all_data), function(i){
  print(paste0("cca association ", i));
  try(retrieve_association_cca(cca_res = cca_all[[i]]$cca_res,
                               index_association_X = all_data[[i]]$index.associated.species,
                               index_association_Y = all_data[[i]]$index.associated.metabolites))
}, mc.cores = 10)
names(cca_conf_tab_all) <- names(all_data)

metrics_cca <- imap(cca_conf_tab_all, ~
                      try(as.data.frame(.x$metrics) %>%
                            mutate(Rep = .y) %>%
                            mutate(Data = Rep %>% str_remove("_.*")) %>%
                            mutate(Method = "sCCA") %>%
                            mutate(Log_metabo = "nolog"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x) %>%
  rowwise() %>% mutate(harmonic.F1 = psych::harmonic.mean(c(Sensitivity, Specificity))) %>%
  pivot_longer(names_to = "Metrics", values_to = "Value", -c(Rep, Data, Method, Log_metabo))

save(metrics_cca, file = "metrics_cca.rda")

ggplot(metrics_cca, aes(x = Method, y = Value)) + geom_violin(aes(fill = Metrics), color = "black") +
  facet_grid(Metrics ~ Data, scales = "free_x") + theme_bw()



# do log for autism and adenomas
all_data_wihout_adenomas <- c(data_autism, data_konzo)

cca_all_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas),function(i){
  print(paste0("scca: ", i));
  OTU = all_data_wihout_adenomas[[i]]$Simulated.Microbiotes;
  METABO = all_data_wihout_adenomas[[i]]$Simulated.Metabolites;

  produce_cca(OTU = OTU,
              METABO = METABO, dolog = TRUE)
}, mc.cores = 50)
names(cca_all_log) <-  names(all_data_wihout_adenomas)


cca_conf_tab_all_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas), function(i){
  print(paste0("cca association ", i));
  try(retrieve_association_cca(cca_res = cca_all_log[[i]]$cca_res,
                               index_association_X = all_data_wihout_adenomas[[i]]$index.associated.species,
                               index_association_Y = all_data_wihout_adenomas[[i]]$index.associated.metabolites))
}, mc.cores = 10)
names(cca_conf_tab_all_log) <- names(all_data_wihout_adenomas)

metrics_cca_log <- imap(cca_conf_tab_all_log, ~
                          try(as.data.frame(.x$metrics) %>%
                                mutate(Rep = .y) %>%
                                mutate(Data = Rep %>% str_remove("_.*")) %>%
                                mutate(Method = "sCCA") %>%
                                mutate(Log_metabo = "log"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x) %>%
  rowwise() %>% mutate(harmonic.F1 = psych::harmonic.mean(c(Sensitivity, Specificity))) %>%
  pivot_longer(names_to = "Metrics", values_to = "Value", -c(Rep, Data, Method, Log_metabo))


ggplot(bind_rows(metrics_cca_log, metrics_cca), aes(x = Method, y = Value)) + geom_violin(aes(fill = Log_metabo), color = "black") +
  facet_grid(Metrics ~ Data, scales = "free_x") + theme_bw()





# penalty param ---------------------------------------------------



bestpenaltyx_nolog <- lapply(cca_all, function(x){
  # autism; adenomas; cca
  lapply(x, function(y)y$perm.out$bestpenaltyx) %>% unlist
}) %>% imap_dfr(~tibble(data = .y, value = .x, mode = "bestpenaltyx", log = FALSE))

bestpenaltyz_nolog <- lapply(cca_all, function(x){
  # autism; adenomas; cca
  lapply(x, function(y)y$perm.out$bestpenaltyz) %>% unlist
}) %>% imap_dfr(~tibble(data = .y, value = .x, mode = "bestpenaltyz",  log = FALSE))


bestpenaltyx_log <- imap_dfr(cca_all_log,
                             ~tibble(data = .y %>% str_remove("_.*"), value = .x$cca$penaltyx, mode = "bestpenaltyx", log = TRUE))

bestpenaltyz_log <- imap_dfr(cca_all_log,
                             ~tibble(data = .y %>% str_remove("_.*"), value = .x$cca$penaltyz, mode = "bestpenaltyz", log = TRUE))

p1 <- rbind(bestpenaltyx_nolog, bestpenaltyz_nolog) %>%
  ggplot(aes(x = value, fill = data)) +
  geom_histogram() +
  facet_grid(data~mode) + theme_bw() + ggtitle(label = "CCA penalty, log = FALSE",subtitle = "penaltyxs=seq(0.01, 1,0.1), penaltyzs=seq(0.01, 1,0.1)")

p2 <- rbind(bestpenaltyx_log, bestpenaltyz_log) %>%
  ggplot(aes(x = value, fill = data)) +
  geom_histogram() +
  facet_grid(data~mode) + theme_bw() + ggtitle(label = "CCA penalty, log = TRUE",subtitle = "penaltyxs=seq(0.01, 1,0.1), penaltyzs=seq(0.01, 1,0.1)")

png("CCA_penalty.png")
p1 + p2
dev.off()
