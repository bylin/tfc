library(glmnet)
library(stringr)
library(ROCR)
library(plyr)

setupTrainingTestingSets <- function(metadata, covariate, counts) {
  # metadata: filtered metadata for covariate of interest (e.g. no [Not Applicable] entries)
  # covariate: name of covariate (must be a column in metadata)
  # counts must be filtered for covariate of interest
  random_sample <- unlist(sapply(levels(metadata[, covariate]), function(class) {
    num_samples <- round(sum(metadata[, covariate] == class)/2)
    sample(metadata[metadata[, covariate] == class, ]$participant_id, num_samples)
  }))
  training_counts <-  counts[, metadata$participant_id %in% random_sample]
  training_metadata <- metadata[metadata$participant_id %in% random_sample, ]
  testing_counts <- counts[, !(metadata$participant_id %in% random_sample)]
  testing_metadata <- metadata[!(metadata$participant_id %in% random_sample), ]
  list(training_counts = training_counts, training_metadata = training_metadata, testing_counts = testing_counts, testing_metadata = testing_metadata, participant_ids = random_sample)
}

buildTestGlm <- function(dataset, features, covariate, randomize = FALSE) {
  training_counts <- t(dataset$training_counts[rownames(dataset$training_counts) %in% features, ])
  training_response <- dataset$training_metadata[match(colnames(dataset$training_counts), dataset$training_metadata$barcode), ][, covariate]
  testing_counts <- t(dataset$testing_counts[rownames(dataset$testing_counts) %in% features, ])
  testing_response <- dataset$testing_metadata[match(colnames(dataset$testing_counts), dataset$testing_metadata$barcode), ][, covariate]
  model <- glmnet(training_counts, training_response, family = "binomial")
  chosen_model <- which.min(abs(model$df - 20)) # select for ~ 20 features)
  pred <- predict(model, newx = testing_counts, s = model$lambda[chosen_model])
  prob <- prediction(pred, testing_response)
  if (randomize) prob <- prediction(pred, sample(testing_response))
  perf <- performance(prob, "tpr", "fpr")
  roc <- data.frame(TPR = unlist(perf@y.values), FPR = unlist(perf@x.values))
  coef <- model$beta[, chosen_model]
  coef <- coef[coef != 0]
  auc <- unlist(performance(prob, "auc")@y.values)
  auc <- round(auc, 3)
  list(roc = roc, coef = coef, auc = auc)
}

parseGlms <- function(glms, class) ldply(1:length(glms), function(i) data.frame(glms[[i]]$roc, Trial = i, AUC = glms[[i]]$auc, Class = class))

plotGlms <- function(roc_df) ggplot(roc_df) + geom_line(aes(x = FPR, y = TPR, color = Class, group = paste0(Trial, Class)), alpha = 0.3) + geom_abline(intercept = 0, slope = 1, colour = "gray") + ylab("TPR") + xlab("FPR") + facet_wrap(~ Class, nrow = 2) + guides(color = FALSE)

plotMeanAucs <- function(all_feature_glms) {
  mean_aucs <- lapply(all_feature_glms, function(glms) {
    aucs <- lapply(1:length(glms), function(i) glms[[i]]$auc)
    mean(unlist(aucs))
  })
  df <- data.frame(AUC = unlist(mean_aucs), Class = as.factor(c("tsRNA", "miRNA", "snoRNA", "piRNA", "tRNA half", "tRNA")))
  ggplot(df) + geom_bar(aes(x = Class, y = AUC), stat = "identity") + ggtitle("Mean AUC by sRNA type")
}
