library(glmnet)
library(stringr)
library(ROCR)
library(plyr)

setupTrainingTestingSets <- function(metadata, covariate, counts) {
  # metadata: filtered metadata for covariate of interest (e.g. no [Not Applicable] entries)
  # covariate: name of covariate (must be a column in metadata)
  # counts must be filtered for covariate of interest
  num_samples_per_class <- min(table(metadata[, covariate]))/2
  random_sample <- unlist(sapply(levels(metadata[, covariate]), function(class) {
    sample(metadata[metadata[, covariate] == class, ]$barcode, num_samples_per_class)
  }))
  training_counts <-  counts[, metadata$barcode %in% random_sample]
  training_metadata <- metadata[metadata$barcode %in% random_sample, ]
  #random_sample_2 <- unlist(sapply(levels(metadata[, covariate]), function(class) {
  #  sample(metadata[metadata[, covariate] == class & !(metadata$barcode %in% random_sample), ]$barcode, num_samples_per_class)
  #}))
  #testing_counts <- counts[, metadata$barcode %in% random_sample_2]
  #testing_metadata <- metadata[metadata$barcode %in% random_sample_2, ]
  testing_counts <- counts[, !(metadata$barcode %in% random_sample)]
  testing_metadata <- metadata[!(metadata$barcode %in% random_sample), ]
  list(training_counts = training_counts, training_metadata = training_metadata, testing_counts = testing_counts, testing_metadata = testing_metadata)
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

extractCoefs <- function(glms, ntrials) {
  df <- ldply(1:ntrials, function(trial) data.frame(name = names(glms[[trial]]$coef), score = 1/length(glms[[trial]]$coef), stringsAsFactors = FALSE))
  df <- dcast(df, name ~ ., fun.aggregate = sum)
  colnames(df) <- c("name", "score")
  df[order(df$score, decreasing = TRUE), ]
}

parseGlms <- function(glms, class) ldply(1:length(glms), function(i) data.frame(glms[[i]]$roc, Trial = i, AUC = glms[[i]]$auc, Class = class))

plotGlms <- function(roc_df) ggplot(roc_df) + geom_line(aes(x = FPR, y = TPR, color = Class, group = paste0(Trial, Class)), alpha = 0.3) + geom_abline(intercept = 0, slope = 1, colour = "gray") + ylab("TPR") + xlab("FPR") + facet_wrap(~ Class, nrow = 2) + guides(color = FALSE)

plotAucs <- function(all_feature_glms) {
  aucs <- ldply(all_feature_glms, function(glms) unlist(lapply(1:length(glms), function(i) glms[[i]]$auc)))
  rownames(aucs) <- c("tsRNA", "miRNA")
  df <- melt(t(aucs))[, -1]
  colnames(df) <- c("sRNA", "AUC")
  ggplot(df) + geom_boxplot(aes(x = sRNA, y = AUC)) + ggtitle("AUC distribution by sRNA type") + ylim(0, 1)
}

buildMultinomialGlm <- function(metadata, counts, covariate, features, randomize = FALSE) {
  # set up training/testing sets. first, set parameters
  dataset <- setupTrainingTestingSets(metadata, covariate, counts)
  training_counts <- t(dataset$training_counts[rownames(dataset$training_counts) %in% features, ])
  training_response <- dataset$training_metadata[match(colnames(dataset$training_counts), dataset$training_metadata$barcode), ][, covariate]
  testing_counts <- t(dataset$testing_counts[rownames(dataset$testing_counts) %in% features, ])
  testing_response <- dataset$testing_metadata[match(colnames(dataset$testing_counts), dataset$testing_metadata$barcode), ][, covariate]
  if (randomize) training_response <- sample(training_response)

  # build and test model
  model <- cv.glmnet(training_counts, training_response, family = "multinomial", type.measure = "class", type.multinomial = "grouped")
  pred <- predict(model, newx = testing_counts, s = "lambda.min", type = "class") 
  unname(confusionMatrix(pred, testing_response)$overall[1])
}

buildCoxGlm <- function(metadata, counts, covariate, features, randomize = FALSE) {
  # set up training/testing sets. first, set parameters
  dataset <- setupTrainingTestingSets(metadata, covariate, counts)
  training_counts <- t(dataset$training_counts[rownames(dataset$training_counts) %in% features, ])
  training_response <- dataset$training_metadata[match(colnames(dataset$training_counts), dataset$training_metadata$barcode), ][, c("days_survived", covariate)]
  testing_counts <- t(dataset$testing_counts[rownames(dataset$testing_counts) %in% features, ])
  testing_response <- dataset$testing_metadata[match(colnames(dataset$testing_counts), dataset$testing_metadata$barcode), ][, covariate]
  if (randomize) training_response <- sample(training_response)

  # build and test model
  model <- cv.glmnet(training_counts, training_response, family = "multinomial", type.measure = "class", type.multinomial = "grouped")
  pred <- predict(model, newx = testing_counts, s = "lambda.min", type = "class") 
  unname(confusionMatrix(pred, testing_response)$overall[1])
}
