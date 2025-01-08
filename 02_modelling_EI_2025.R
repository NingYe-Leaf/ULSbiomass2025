library(readxl)
library(missForest)
library(caret)
library(randomForest)
library(ggplot2)
library(xgboost)
library(corrplot)
library(gridExtra)
library(openxlsx)
library(treeshap)
library(shapviz)
library(pls)

######################################################################################
data <- read_excel("G:/Project2023/03_Analysis/combined_data_0.1_AGB.xls")
data_subset <- data[, c(4, 41:137)]#AGB
colnames(data_subset)
data_subset$Species <- as.factor(data_subset$Species)
dummy_species <- model.matrix(~ Species - 1, data = data_subset)
data_subset <- cbind(data_subset, dummy_species)
data_subset <- data_subset[, c(2:100)]
colnames(data_subset)
data_subset <- as.data.frame(lapply(data_subset, as.numeric))

numeric_columns <- sapply(data_subset, is.numeric)
data_subset[, numeric_columns] <- lapply(data_subset[, numeric_columns], function(x) ifelse(!is.finite(x), NA, x))
threshold <- nrow(data_subset) * 0.1

removed_columns <- character(0)
for (col in names(data_subset)) {
  na_count <- sum(is.na(data_subset[[col]]))
  if (na_count >= threshold) {
    data_subset[[col]] <- NULL
    removed_columns <- c(removed_columns, col)
    cat(paste("Removed column:", col, "\n"))
  }
}

cat("\nColumn names with less than 10% NA values:\n")
selected_columns <- names(data_subset)[colSums(is.na(data_subset)) < threshold & colSums(is.na(data_subset)) > 0]
print(selected_columns)

set.seed(123)  
data_subset<- missForest(data_subset)
imputed_data <- data_subset$ximp

std_dev <- apply(imputed_data, 2, function(x) sd(x, na.rm = TRUE))
selected_columns <- imputed_data[, !is.na(std_dev) & std_dev > 0]
colnames(selected_columns)
selected_columnsnew <- selected_columns[, c(1:86)]
colnames(selected_columnsnew)
correlation_matrix <- cor(selected_columnsnew)

corrplot(correlation_matrix, method = 'circle', type = 'upper', insig = 'blank',
         order = 'AOE', diag = FALSE,  tl.cex = 0.6, tl.col = "black")$corrPos -> p1

set.seed(123) 
train_indices <- sample(1:nrow(selected_columns), 0.6 * nrow(selected_columns))
train_data <- selected_columns[train_indices, ]
test_data <- selected_columns[-train_indices, ]

write.xlsx(train_data, "G:/Project2023/03_Analysis/traindata_AGB.xlsx", sheetName = "Sheet1", tableStyle = "TableStyleMedium9", colNames = TRUE)

############################# random forest#######################################
customRF <- list(type = "Regression",
                 library = "randomForest",
                 loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree","nodesize"),
                                  class = rep("numeric", 3),
                                  label = c("mtry", "ntree", "Nodesize"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree = param$ntree,
               nodesize = param$nodesize,#set node size
               importance = TRUE)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "response")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

control <- trainControl(method = "cv", number = 10)
tunegrid <- expand.grid(mtry = 1:44, ntree = c(500, 1000, 1500, 2000, 2500),nodesize = c(5, 10, 15, 20))

set.seed(123)
start <- Sys.time()
custom <- train(AGB ~ ., data = train_data, 
                method = customRF, 
                metric = "RMSE", 
                tuneGrid = tunegrid,
                trControl = control)
end <- Sys.time()
end-start

summary(custom)
plot(custom)

best_custom_model <- custom$finalModel
custom$bestTune
print(best_custom_model)
plot(best_custom_model)
plot(best_custom_model$rsq)
plot(best_custom_model$mse)
#The varImp() function comes from the caret package.
#The varImp() function in caret attempts to offer a more standardized way to view variable importance across various models.
var_importance <- varImp(best_custom_model)
print(var_importance)
#The importance() function originates from the randomForest package.
#With the importance() function, the return is a matrix or data frame that contains quantitative measures of importance for each 
#variable in the random forest model. These measures may include Mean Decrease Accuracy and/or Mean Decrease Gini, among others.
var_importance2 <-importance(best_custom_model)
print(var_importance2)

#####################################
######################################
best_custom_model_shap<-randomForest.unify(best_custom_model, train_data)
treeshap_res <- treeshap(best_custom_model_shap, train_data)
treeshap_res$shaps
plot_contribution(treeshap_res,max_vars = 10,obs = 1)
plot_feature_importance(treeshap_res,desc_sorting = TRUE,max_vars = 10,title = "Feature Importance",subtitle = NULL)
######################################
######################################

# Create a variable importance plot
ggplot(importance_data, aes(x = reorder(rownames(importance_data), -Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "blue") +
  xlab("Variables") +
  ylab("Variable Importance") +
  ggtitle("Variable Importance Plot") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate R-squared and RMSE for test_data_plot
test_data_plot<-test_data
test_data_plot$predicted_AGB <- predict(best_custom_model, newdata = test_data_plot)
test_data_plot$residuals <- test_data_plot$AGB - test_data_plot$predicted_AGB
ss_total_test <- sum((test_data_plot$AGB - mean(test_data_plot$AGB))^2)
ss_residual_test <- sum(test_data_plot$residuals^2)
r_squared_test <- 1 - (ss_residual_test / ss_total_test)
rmse_test <- sqrt(mean(test_data_plot$residuals^2))

cat("RMSE for test data:", rmse_test, "\n")
cat("Rsquared for test data:", r_squared_test, "\n")

# Create a new data frame with selected columns and the filter condition
globoidea_data <- subset(test_data_plot, SpeciesE..globoidea == 1, select = c(AGB, predicted_AGB, SpeciesE..globoidea))
globoidea_data$residuals <- globoidea_data$AGB - globoidea_data$predicted_AGB
ss_total_globoidea_data <- sum((globoidea_data$AGB - mean(globoidea_data$AGB))^2)
ss_residual_globoidea_data <- sum(globoidea_data$residuals^2)
r_squared_globoidea_data <- 1 - (ss_residual_globoidea_data / ss_total_globoidea_data)
rmse_globoidea_data <- sqrt(mean(globoidea_data$residuals^2))

cat("RMSE for globoidea:", rmse_globoidea_data, "\n")
cat("Rsquared for globoidea:", r_squared_globoidea_data, "\n")

bosistoana_data <- subset(test_data_plot, SpeciesE..bosistoana == 1, select = c(AGB, predicted_AGB, SpeciesE..bosistoana))
bosistoana_data$residuals <- bosistoana_data$AGB - bosistoana_data$predicted_AGB
ss_total_bosistoana_data <- sum((bosistoana_data$AGB - mean(bosistoana_data$AGB))^2)
ss_residual_bosistoana_data <- sum(bosistoana_data$residuals^2)
r_squared_bosistoana_data <- 1 - (ss_residual_bosistoana_data / ss_total_bosistoana_data)
rmse_bosistoana_data <- sqrt(mean(bosistoana_data$residuals^2))

cat("RMSE for bosistoana:", rmse_bosistoana_data, "\n")
cat("Rsquared for bosistoana:", r_squared_bosistoana_data, "\n")

# Create plot for test_data
plot_test_data <- ggplot(test_data_plot, aes(x = AGB, y = predicted_AGB)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual AGB (kg)", y = "Predicted AGB (kg)") +
  ggtitle("Test Data Predicted vs. Actual AGB") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

# Create a residual plot for test_data
residual_plot_test <- ggplot(test_data_plot, aes(x = predicted_AGB, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Actual AGB (kg)", y = "Residuals (kg)") +
  ggtitle("All test set (41 trees)") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

grid.arrange(plot_test_data, residual_plot_test, ncol = 2)  

# Create a new column to determine point color based on conditions
test_data_plot$point_color <- ifelse(test_data_plot$SpeciesE..bosistoana == 1, "SpeciesE..bosistoana",
                                     ifelse(test_data_plot$SpeciesE..globoidea == 1, "SpeciesE..globoidea", "Other"))

# Create plot for test_data with colored points
plot_test_data <- ggplot(test_data_plot, aes(x = AGB, y = predicted_AGB, color = point_color)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("SpeciesE..bosistoana" = "blue", "SpeciesE..globoidea" = "orange", "Other" = "black")) +
  labs(x = "Actual AGB (kg)", y = "Predicted AGB (kg)") +
  ggtitle("Test Data Predicted vs. Actual AGB") +
  theme(panel.background = element_rect(fill = "white", color = "black"))
print(plot_test_data)

# Create one residual plot for different species
residual_plot <- ggplot(test_data_plot, aes(x = test_data_plot$AGB, y = test_data_plot$residuals, color = point_color)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("SpeciesE..bosistoana" = "blue", "SpeciesE..globoidea" = "orange", "Other" = "black")) +
  labs(x = "Actual AGB (kg)", y = "Residuals") +
  ggtitle("Residual Plot by Species") +
  theme(panel.background = element_rect(fill = "white", color = "black"))
print(residual_plot)
 
# Create two residual plot for different species
residual_plot_bosistoana <- ggplot(bosistoana_data, aes(x = predicted_AGB, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted AGB (kg)", y = "Residuals (kg)") +
  ggtitle("E. bosistoana in test set (32 trees)") +
  theme(panel.background = element_rect(fill = "white", color = "black"))


residual_plot_globoidea <- ggplot(globoidea_data, aes(x = predicted_AGB, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted AGB (kg)", y = "Residuals (kg)") +
  ggtitle("E. globoidea in test set (9 trees)") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

# Display the plots and residual plots side by side
grid.arrange( residual_plot_bosistoana, residual_plot_globoidea , ncol = 2)

##############################XGBoost###################################
param_grid <- expand.grid(
  nrounds = c(50, 100, 150),
  max_depth = c(2, 4, 6, 8, 10),  # Maximum depth of trees
  eta = c(0.01, 0.05, 0.1, 0.2, 0.3),  # Fraction of features used in each boosting round
  min_child_weight = c(1, 2, 3, 4, 5),  # Minimum sum of instance weight (hessian) needed in a child
  subsample = c(0.6, 0.7, 0.8, 0.9, 1),  # Fraction of training data to randomly sample in each boosting round
  gamma =c(0, 0.1, 0.2, 0.3, 0.4)
)
dtrain <- xgb.DMatrix(data.matrix(train_data[, -which(names(train_data) == "AGB")]), label = train_data$AGB)
start <- Sys.time()
results <- list()
set.seed(123)
for (i in 1:nrow(param_grid)) {
  params <- list(
    booster = "gbtree",
    objective = "reg:squarederror",
    nthread = 2,
    eval_metric = "rmse",
    max_depth = param_grid$max_depth[i],
    eta = param_grid$eta[i],
    min_child_weight = param_grid$min_child_weight[i],
    subsample = param_grid$subsample[i],
    gamma = param_grid$gamma[i]
  )
  cv <- xgb.cv(nrounds = param_grid$nrounds[i], params = params, data = dtrain, nfold = 10)
  results[[i]] <- cv$evaluation_log
}

best_model_index <- which.min(sapply(results, function(log) min(log$test_rmse_mean)))
best_params <- param_grid[best_model_index, ]
cat("Best parameters:")
print(best_params)

end <- Sys.time()
end-start

final_params <- list(
  booster = "gbtree",
  objective = "reg:squarederror",
  nthread = 2,
  eval_metric = "rmse",
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  min_child_weight = best_params$min_child_weight,
  subsample = best_params$subsample,
  gamma = best_params$gamma
)

set.seed(123)
best_xgb_model <- xgb.train(
  params = final_params,
  data = dtrain,
  nrounds = best_params$nrounds,  
  verbose = 1  
)

####################shapely############################################
xgb.importance(colnames(dtrain), model = best_xgb_model)
shp <- shapviz(best_xgb_model, X_pred = data.matrix(train_data[, -which(names(train_data) == "AGB")]), X = train_data)
sv_waterfall(shp, row_id = 1,max_display = 11L)
sv_importance(shp)
sv_importance(shp, kind = "beeswarm", alpha = 0.2, width = 0.2,max_display = 10L)
sv_importance(shp, kind = "both", alpha = 0.4, width = 0.2,max_display = 10L)
sv_importance(shp, kind = "beeswarm",max_display = 10L,bee_width = 0.5)

best_xgb_model_shap<-xgboost.unify(best_xgb_model, train_data)
treeshap_res <- treeshap(best_xgb_model_shap, train_data,interactions = TRUE)
treeshap_res$shaps
plot_contribution(treeshap_res,max_vars = 10)
plot_feature_dependence(treeshap_res,variable="HOME",title = "Feature Dependence",subtitle = NULL)
plot_feature_importance(treeshap_res,desc_sorting = TRUE,max_vars = 10,title = "Feature Importance",subtitle = NULL)
plot_interaction(treeshap_res,var1="zmax",var2="vzrump",title = "SHAP Interaction Value Plot",subtitle =  NULL)


######################################
dtest <- xgb.DMatrix(data.matrix(test_data[, -which(names(test_data) == "AGB")]))
xgb_predictions <- predict(best_xgb_model, dtest)

xgb_actual_values <- test_data$AGB
xgb_rmse <- sqrt(mean((xgb_predictions - xgb_actual_values)^2))

xgb_ss_total <- sum((xgb_actual_values - mean(xgb_actual_values))^2)
xgb_ss_residual <- sum((xgb_actual_values - xgb_predictions)^2)
xgb_r_squared <- 1 - (xgb_ss_residual / xgb_ss_total)

cat("XGBoost RMSE:", xgb_rmse, "\n")
cat("XGBoost R-squared:", xgb_r_squared, "\n")

test_data$predicted_AGB <- xgb_predictions
test_data$residuals <- test_data$AGB - test_data$predicted_AGB

globoidea_data <- subset(test_data, SpeciesE..globoidea == 1, select = c(AGB, predicted_AGB, SpeciesE..globoidea))
globoidea_data$residuals <- globoidea_data$AGB - globoidea_data$predicted_AGB
ss_total_globoidea_data <- sum((globoidea_data$AGB - mean(globoidea_data$AGB))^2)
ss_residual_globoidea_data <- sum(globoidea_data$residuals^2)
r_squared_globoidea_data <- 1 - (ss_residual_globoidea_data / ss_total_globoidea_data)
rmse_globoidea_data <- sqrt(mean(globoidea_data$residuals^2))

cat("RMSE for globoidea:", rmse_globoidea_data, "\n")
cat("Rsquared for globoidea:", r_squared_globoidea_data, "\n")

bosistoana_data <- subset(test_data, SpeciesE..bosistoana == 1, select = c(AGB, predicted_AGB, SpeciesE..bosistoana))
bosistoana_data$residuals <- bosistoana_data$AGB - bosistoana_data$predicted_AGB
ss_total_bosistoana_data <- sum((bosistoana_data$AGB - mean(bosistoana_data$AGB))^2)
ss_residual_bosistoana_data <- sum(bosistoana_data$residuals^2)
r_squared_bosistoana_data <- 1 - (ss_residual_bosistoana_data / ss_total_bosistoana_data)
rmse_bosistoana_data <- sqrt(mean(bosistoana_data$residuals^2))

cat("RMSE for bosistoana:", rmse_bosistoana_data, "\n")
cat("Rsquared for bosistoana:", r_squared_bosistoana_data, "\n")

plot_test_data <- ggplot(test_data, aes(x = AGB, y = predicted_AGB)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual AGB (kg)", y = "Predicted AGB (kg)") +
  ggtitle("Test Data Predicted vs. Actual AGB") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

residual_plot_test <- ggplot(test_data, aes(x = predicted_AGB, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted AGB (kg)", y = "Residuals (kg)") +
  ggtitle("All test set (41 trees)") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

grid.arrange(plot_test_data, residual_plot_test, ncol = 2) 

test_data$point_color <- ifelse(test_data$SpeciesE..bosistoana == 1, "SpeciesE..bosistoana",
                                     ifelse(test_data$SpeciesE..globoidea == 1, "SpeciesE..globoidea", "Other"))
plot_test_data <- ggplot(test_data, aes(x = AGB, y = predicted_AGB, color = point_color)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("SpeciesE..bosistoana" = "blue", "SpeciesE..globoidea" = "orange", "Other" = "black")) +
  labs(x = "Actual AGB (kg)", y = "Predicted AGB (kg)") +
  ggtitle("Test Data Predicted vs. Actual AGB") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

print(plot_test_data)

residual_plot <- ggplot(test_data, aes(x = test_data$AGB, y = test_data$residuals, color = point_color)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("SpeciesE..bosistoana" = "blue", "SpeciesE..globoidea" = "orange", "Other" = "black")) +
  labs(x = "Actual AGB (kg)", y = "Residuals") +
  ggtitle("Residual Plot by Species") +
  theme(panel.background = element_rect(fill = "white", color = "black"))
print(residual_plot)

residual_plot_bosistoana <- ggplot(bosistoana_data, aes(x = predicted_AGB, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted AGB (kg)", y = "Residuals (kg)") +
  ggtitle("E. bosistoana in test set (32 trees)") +
  theme(panel.background = element_rect(fill = "white", color = "black"))


residual_plot_globoidea <- ggplot(globoidea_data, aes(x = predicted_AGB, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted AGB (kg)", y = "Residuals (kg)") +
  ggtitle("E. globoidea in test set (9 trees)") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

grid.arrange(residual_plot_bosistoana, residual_plot_globoidea , ncol = 2)


####################PLSR############################################
start <- Sys.time()
set.seed(123)
ctrl <- trainControl(method = "cv", number = 10)
plsTune <- train(AGB ~ ., data = train_data,
                 method = "pls", 
                 # set hyper-parameter tuning range
                 tuneGrid = expand.grid(.ncomp = 1:50),
                 trControl = ctrl)
end <- Sys.time()
end-start
print(plsTune)
plot(plsTune, main = 'PLSR Model Tuning: Number of Components vs Performance')
plsImp <- varImp(plsTune, scale = FALSE)
plot(plsImp, top = 10, scales = list(y = list(cex = 0.95)))

final_ncomp <- plsTune$bestTune$ncomp
final_ncomp
set.seed(123)
final_plsr_model <- plsr(AGB ~ ., data = train_data, ncomp = final_ncomp, scale = FALSE)
print(final_plsr_model)
 
test_data_plsr<-test_data[,1:89]
test_data_plsr$predicted_AGB <- predict(final_plsr_model , newdata = test_data_plsr,ncomp = final_ncomp)

performance_metrics <- postResample(pred = test_data_plsr$predicted_AGB, obs = test_data_plsr$AGB)
cat("RMSE for test data:", performance_metrics['RMSE'], "\n")
cat("Rsquared for test data:", performance_metrics['Rsquared'], "\n")

test_data_plsr$residuals <- test_data_plsr$AGB - test_data_plsr$predicted_AGB
test_data_plsr$point_color <- ifelse(test_data_plsr$SpeciesE..bosistoana == 1, "SpeciesE..bosistoana",
                                ifelse(test_data_plsr$SpeciesE..globoidea == 1, "SpeciesE..globoidea", "Other"))

plot_test_data <- ggplot(test_data_plsr, aes(x = AGB, y = predicted_AGB, color = point_color)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("SpeciesE..bosistoana" = "blue", "SpeciesE..globoidea" = "orange", "Other" = "black")) +
  labs(x = "Actual AGB (kg)", y = "Predicted AGB (kg)") +
  ggtitle("Test Data Predicted vs. Actual AGB") +
  theme(panel.background = element_rect(fill = "white", color = "black"))

print(plot_test_data)

residual_plot <- ggplot(test_data_plsr, aes(x = test_data_plsr$AGB, y = test_data_plsr$residuals, color = point_color)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("SpeciesE..bosistoana" = "blue", "SpeciesE..globoidea" = "orange", "Other" = "black")) +
  labs(x = "Actual AGB (kg)", y = "Residuals") +
  ggtitle("Residual Plot by Species") +
  theme(panel.background = element_rect(fill = "white", color = "black"))
print(residual_plot)

residuals_globoidea <- subset(test_data_plsr, SpeciesE..globoidea == 1)$residuals
actual_values_globoidea <- subset(test_data_plsr, SpeciesE..globoidea == 1)$AGB

residuals_bosistoana <- subset(test_data_plsr, SpeciesE..bosistoana == 1)$residuals
actual_values_bosistoana <- subset(test_data_plsr, SpeciesE..bosistoana == 1)$AGB

rmse_globoidea <- sqrt(mean(residuals_globoidea^2))
ss_total_globoidea <- sum((actual_values_globoidea - mean(actual_values_globoidea))^2)
ss_residual_globoidea <- sum(residuals_globoidea^2)
r_squared_globoidea <- 1 - (ss_residual_globoidea / ss_total_globoidea)

rmse_bosistoana <- sqrt(mean(residuals_bosistoana^2))
ss_total_bosistoana <- sum((actual_values_bosistoana - mean(actual_values_bosistoana))^2)
ss_residual_bosistoana <- sum(residuals_bosistoana^2)
r_squared_bosistoana <- 1 - (ss_residual_bosistoana / ss_total_bosistoana)

# Print the results
cat("RMSE for globoidea:", rmse_globoidea, "\n")
cat("Rsquared for globoidea:", r_squared_globoidea, "\n")

cat("RMSE for bosistoana:", rmse_bosistoana, "\n")
cat("Rsquared for bosistoana:", r_squared_bosistoana, "\n")
