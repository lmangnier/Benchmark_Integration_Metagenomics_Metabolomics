library(mixOmics)
library(compositions)
library(tidyverse)

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
  nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
  if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
  return(x)
}

apply_tune_spls <- function(OTU, METABO, mode = "canonical", list.keepX = NULL,
                            list.keepY = NULL, regression_mode = NULL,
                            dolog = FALSE
){
  # CLR OTU
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

  if(purrr::is_null(list.keepX)){
    list.keepX <- round(seq(from = ncol(X)*0.02, to = ncol(X)*0.1, length.out = 5))
  }
  if(purrr::is_null(list.keepY)){
    list.keepY <- round(seq(from = ncol(Y)*0.02, to = ncol(Y)*0.1, length.out = 5))
  }

  final.spls <- list()
  if(mode == "canonical"){
    set.seed(123)
    tune_res <- try(mixOmics::tune.spls(X, Y, ncomp = 2,
                                        test.keepX = list.keepX,
                                        test.keepY = list.keepY,
                                        nrepeat = 1, folds = 10, # use 10 folds
                                        mode = 'canonical', measure="cor")
    )
    if(is(tune_res) != "try-error"){
      # overwrite final.spls or keep empty
      final.spls = spls(X,  Y,
                        ncomp = tune_res$choice.ncomp,
                        keepX = tune_res$choice.keepX[tune_res$choice.ncomp],
                        keepY = tune_res$choice.keepY[tune_res$choice.ncomp],
                        mode = "canonical")
    }
  } else if(mode == "regression" & regression_mode == 1){
    # regression X = X, Y = Y
    set.seed(123)
    tune_res <- try(mixOmics::tune.spls(X = X, Y = Y, ncomp = 2,
                                        test.keepX = list.keepX,
                                        test.keepY = list.keepY,
                                        nrepeat = 1, folds = 10, # use 10 folds
                                        mode = 'regression', measure="cor"))
    if(is(tune_res) != "try-error"){
      final.spls = spls(X = X,  Y = Y,
                        ncomp = tune_res$choice.ncomp,
                        keepX = tune_res$choice.keepX[tune_res$choice.ncomp],
                        keepY = tune_res$choice.keepY[tune_res$choice.ncomp],
                        mode = "regression")}
  } else if(mode == "regression" & regression_mode == 2){
    # regression 2 X = Y, Y = X
    set.seed(123)
    tune_res <- try(mixOmics::tune.spls(X = Y, Y = X, ncomp = 2,
                                        test.keepX = list.keepY, # swipe keepX and keepX
                                        test.keepY = list.keepX,
                                        nrepeat = 1, folds = 10, # use 10 folds
                                        mode = 'regression', measure="cor"))
    if(is(tune_res) != "try-error"){
      final.spls  <- spls(X = Y,  Y = X,
                          ncomp = tune_res$choice.ncomp,
                          keepX = tune_res$choice.keepX[tune_res$choice.ncomp],
                          keepY = tune_res$choice.keepY[tune_res$choice.ncomp],
                          mode = "regression")}
  }
  # can be empty
  return(final.spls)
}


retrieve_association <- function(sPLS_res, index_association_X, index_association_Y){

  ncomp <- sPLS_res$ncomp
  selected_features_X <- lapply(1:ncomp, function(comp){
    mixOmics::selectVar(sPLS_res, comp = comp, block = "X")$X$name}) %>%
    unlist()

  selected_features_Y <- lapply(1:ncomp, function(comp){
    mixOmics::selectVar(sPLS_res, comp = comp, block = "Y")$Y$name}) %>%
    unlist()

  colnames_index_association_X <- paste0("X_", index_association_X)
  colnames_index_association_Y <- paste0("Y_", index_association_Y)

  binded <- data.frame(feature = c(sPLS_res$names$colnames$X, sPLS_res$names$colnames$Y)) %>%
    mutate(truth = ifelse(feature %in% c(colnames_index_association_X, colnames_index_association_Y), 1, 0)) %>%
    mutate(pred = ifelse(feature %in% c(selected_features_X, selected_features_Y), 1, 0)) %>%
    mutate(block = str_extract(feature, "^."))

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


data_autism <- readRDS("synthetic_data_Autism_alter_with_index_associations.RDS")
data_adenomas <- readRDS("synthetic_data_Adenomas_alter_with_index_associations.RDS")
data_konzo <- readRDS("synthetic_data_from_Konzo_dataset_alter_with_index_associations_full.RDS")

names(data_autism) <- paste0("Autism_", seq_along(data_autism))
names(data_adenomas) <- paste0("Adenomas_", seq_along(data_adenomas))
names(data_konzo) <- paste0("Konzo_", seq_along(data_konzo))

all_data <- c(data_autism, data_adenomas, data_konzo)

spls_all_canonical <- parallel::mclapply(seq_along(all_data),function(i){
  print(paste0("spls all canonical: ", i));
  OTU = all_data[[i]]$Simulated.Microbiotes;
  METABO = all_data[[i]]$Simulated.Metabolites;
  list.keepX = round(seq(from = ncol(OTU)*0.02, to = ncol(OTU)*0.2, length.out = 5));
  list.keepY = round(seq(from = ncol(METABO)*0.02, to = ncol(METABO)*0.2, length.out = 5));

  apply_tune_spls(OTU = OTU,
                  METABO = METABO,
                  list.keepX = list.keepX,
                  list.keepY = list.keepY, mode = "canonical")
}, mc.cores = 50)

spls_all_regression_1 <- parallel::mclapply(seq_along(all_data),function(i){
  print(paste0("spls all regression 1: ", i));
  OTU = all_data[[i]]$Simulated.Microbiotes;
  METABO = all_data[[i]]$Simulated.Metabolites;
  list.keepX = round(seq(from = ncol(OTU)*0.02, to = ncol(OTU)*0.2, length.out = 5));
  list.keepY = round(seq(from = ncol(METABO)*0.02, to = ncol(METABO)*0.2, length.out = 5));

  apply_tune_spls(OTU = OTU,
                  METABO = METABO,
                  list.keepX = list.keepX,
                  list.keepY = list.keepY, mode = "regression", regression_mode = 1)
}, mc.cores = 50)

spls_all_regression_2 <- parallel::mclapply(seq_along(all_data),function(i){
  print(paste0("spls all regression 2: ", i));
  OTU = all_data[[i]]$Simulated.Microbiotes;
  METABO = all_data[[i]]$Simulated.Metabolites;
  list.keepX = round(seq(from = ncol(OTU)*0.02, to = ncol(OTU)*0.2, length.out = 5));
  list.keepY = round(seq(from = ncol(METABO)*0.02, to = ncol(METABO)*0.2, length.out = 5));

  apply_tune_spls(OTU = OTU,
                  METABO = METABO,
                  list.keepX = list.keepX,
                  list.keepY = list.keepY, mode = "regression", regression_mode = 2)
}, mc.cores = 50)

save(spls_all_canonical, spls_all_regression_1, spls_all_regression_2, file = "spls_all_res.rda")


spls_conf_table_all_canonical <- parallel::mclapply(seq_along(all_data), function(i){
  print(paste0("spls canonical association ", i));
  try(retrieve_association(sPLS_res = spls_all_canonical[[i]],
                           index_association_X = all_data[[i]]$index.associated.species,
                           index_association_Y = all_data[[i]]$index.associated.metabolites))
}, mc.cores = 10)

spls_conf_table_all_regression_1 <- parallel::mclapply(seq_along(all_data), function(i){
  print(paste0("spls regrssion 1 ", i));
  try(retrieve_association(sPLS_res = spls_all_regression_1[[i]],
                           index_association_X = all_data[[i]]$index.associated.species,
                           index_association_Y = all_data[[i]]$index.associated.metabolites))
}, mc.cores = 10)

spls_conf_table_all_regression_2 <- parallel::mclapply(seq_along(all_data), function(i){
  print(paste0("spls regrssion 2 ", i));
  try(retrieve_association(sPLS_res = spls_all_regression_2[[i]],
                           index_association_X = all_data[[i]]$index.associated.metabolites,
                           index_association_Y = all_data[[i]]$index.associated.species))
}, mc.cores = 10)

names(spls_conf_table_all_canonical) <- names(spls_conf_table_all_regression_1) <- names(spls_conf_table_all_regression_2) <- names(all_data)
# save(spls_conf_table_all_canonical, spls_conf_table_all_regression_1, spls_conf_table_all_regression_2, file = "spls_conf_tab_final.rda")
#
# load("spls_conf_tab_final.rda")
metrics_spls_canonical <- imap(spls_conf_table_all_canonical, ~
                                 try(as.data.frame(.x$metrics) %>%
                                       mutate(Rep = .y) %>%
                                       mutate(Data = Rep %>% str_remove("_.*")) %>%
                                       mutate(Method = "sPLS-Can"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x)

metrics_spls_regression_1 <- imap(spls_conf_table_all_regression_1, ~
                                    try(as.data.frame(.x$metrics) %>%
                                          mutate(Rep = .y) %>%
                                          mutate(Data = Rep %>% str_remove("_.*")) %>%
                                          mutate(Method = "sPLS-Reg1"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x)

metrics_spls_regression_2 <- imap(spls_conf_table_all_regression_2, ~
                                    try(as.data.frame(.x$metrics) %>%
                                          mutate(Rep = .y) %>%
                                          mutate(Data = Rep %>% str_remove("_.*")) %>%
                                          mutate(Method = "sPLS-Reg2"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x)

metrics_spls <- bind_rows(metrics_spls_canonical, metrics_spls_regression_1, metrics_spls_regression_2) %>%
  rowwise() %>% mutate(harmonic.F1 = psych::harmonic.mean(c(Sensitivity, Specificity))) %>%
  pivot_longer(names_to = "Metrics", values_to = "Value", -c(Rep, Data, Method))

# save(metrics_spls, file = "metrics_spls.rda")
ggplot(metrics_spls, aes(x = Method, y = Value)) + geom_violin(aes(fill = Metrics), color = "black") +
  facet_grid(Metrics ~ Data, scales = "free_x") + theme_bw()


# do log on
all_data_wihout_adenomas <- c(data_autism, data_konzo)

spls_all_canonical_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas),function(i){
  print(paste0("spls all canonical: ", i));
  OTU = all_data_wihout_adenomas[[i]]$Simulated.Microbiotes;
  METABO = all_data_wihout_adenomas[[i]]$Simulated.Metabolites;
  list.keepX = round(seq(from = ncol(OTU)*0.02, to = ncol(OTU)*0.2, length.out = 5));
  list.keepY = round(seq(from = ncol(METABO)*0.02, to = ncol(METABO)*0.2, length.out = 5));

  apply_tune_spls(OTU = OTU,
                  METABO = METABO,
                  list.keepX = list.keepX,
                  list.keepY = list.keepY, mode = "canonical", dolog = TRUE)
}, mc.cores = 50)

spls_all_regression_1_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas),function(i){
  print(paste0("spls all regression 1: ", i));
  OTU = all_data_wihout_adenomas[[i]]$Simulated.Microbiotes;
  METABO = all_data_wihout_adenomas[[i]]$Simulated.Metabolites;
  list.keepX = round(seq(from = ncol(OTU)*0.02, to = ncol(OTU)*0.2, length.out = 5));
  list.keepY = round(seq(from = ncol(METABO)*0.02, to = ncol(METABO)*0.2, length.out = 5));

  apply_tune_spls(OTU = OTU,
                  METABO = METABO,
                  list.keepX = list.keepX,
                  list.keepY = list.keepY, mode = "regression", regression_mode = 1, dolog = TRUE)
}, mc.cores = 50)

spls_all_regression_2_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas),function(i){
  print(paste0("spls all regression 2: ", i));
  OTU = all_data_wihout_adenomas[[i]]$Simulated.Microbiotes;
  METABO = all_data_wihout_adenomas[[i]]$Simulated.Metabolites;
  list.keepX = round(seq(from = ncol(OTU)*0.02, to = ncol(OTU)*0.2, length.out = 5));
  list.keepY = round(seq(from = ncol(METABO)*0.02, to = ncol(METABO)*0.2, length.out = 5));

  apply_tune_spls(OTU = OTU,
                  METABO = METABO,
                  list.keepX = list.keepX,
                  list.keepY = list.keepY, mode = "regression", regression_mode = 2, dolog = TRUE)
}, mc.cores = 50)

spls_conf_table_all_canonical_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas), function(i){
  print(paste0("spls canonical association ", i));
  try(retrieve_association(sPLS_res = spls_all_canonical_log[[i]],
                           index_association_X = all_data_wihout_adenomas[[i]]$index.associated.species,
                           index_association_Y = all_data_wihout_adenomas[[i]]$index.associated.metabolites))
}, mc.cores = 10)

spls_conf_table_all_regression_1_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas), function(i){
  print(paste0("spls regrssion 1 ", i));
  try(retrieve_association(sPLS_res = spls_all_regression_1_log[[i]],
                           index_association_X = all_data_wihout_adenomas[[i]]$index.associated.species,
                           index_association_Y = all_data_wihout_adenomas[[i]]$index.associated.metabolites))
}, mc.cores = 10)

spls_conf_table_all_regression_2_log <- parallel::mclapply(seq_along(all_data_wihout_adenomas), function(i){
  print(paste0("spls regrssion 2 ", i));
  try(retrieve_association(sPLS_res = spls_all_regression_2_log[[i]],
                           index_association_X = all_data_wihout_adenomas[[i]]$index.associated.metabolites,
                           index_association_Y = all_data_wihout_adenomas[[i]]$index.associated.species))
}, mc.cores = 10)

names(spls_conf_table_all_canonical_log) <- names(spls_conf_table_all_regression_1_log) <- names(spls_conf_table_all_regression_2_log) <- names(all_data_wihout_adenomas)
# save(spls_conf_table_all_canonical_log, spls_conf_table_all_regression_1_log, spls_conf_table_all_regression_2_log, file = "conf_tab_spls_log_all.Rda")
#
# load("conf_tab_spls_log_all.Rda")
metrics_spls_canonical_log <- imap(spls_conf_table_all_canonical_log, ~
                                     try(as.data.frame(.x$metrics) %>%
                                           mutate(Rep = .y) %>%
                                           mutate(Data = Rep %>% str_remove("_.*")) %>%
                                           mutate(Method = "sPLS-Can"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x)

metrics_spls_regression_1_log <- imap(spls_conf_table_all_regression_1_log, ~
                                        try(as.data.frame(.x$metrics) %>%
                                              mutate(Rep = .y) %>%
                                              mutate(Data = Rep %>% str_remove("_.*")) %>%
                                              mutate(Method = "sPLS-Reg1"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x)

metrics_spls_regression_2_log <- imap(spls_conf_table_all_regression_2_log, ~
                                        try(as.data.frame(.x$metrics) %>%
                                              mutate(Rep = .y) %>%
                                              mutate(Data = Rep %>% str_remove("_.*")) %>%
                                              mutate(Method = "sPLS-Reg2"))) %>%
  purrr::discard(~is(.x, "try-error")) %>%
  imap_dfr(~.x)

metrics_spls_log <- bind_rows(metrics_spls_canonical_log, metrics_spls_regression_1_log, metrics_spls_regression_2_log) %>%
  rowwise() %>% mutate(harmonic.F1 = psych::harmonic.mean(c(Sensitivity, Specificity))) %>%
  pivot_longer(names_to = "Metrics", values_to = "Value", -c(Rep, Data, Method))

# save(metrics_spls_log, file = "metrics_spls_log.rda")
ggplot(metrics_spls_log, aes(x = Method, y = Value)) + geom_violin(aes(fill = Metrics), color = "black") +
  facet_grid(Metrics ~ Data, scales = "free_x") + theme_bw() + ggtitle("sPLS, log metabolite")


