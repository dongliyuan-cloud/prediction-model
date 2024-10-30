library(survival)
library(dplyr)
library(writexl)
library(survminer)
library(boot)
library(timeROC)
library(rms)
library(ggplot2)
library(nomogramFormula)

# Setting up a working directory and reading data
setwd("C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2")    ### your working_ directory
data_2 <- read.csv("430.csv", fileEncoding = "GB2312")

#####
str(data_2)
summary(data_2)  

# Converting dichotomous variables to factors and setting reference levels
data_2$sex <- relevel(factor(data_2$sex), ref = "0") 
data_2$smoke <- relevel(factor(data_2$smoke), ref = "0")
data_2$alcohol <- relevel(factor(data_2$alcohol), ref = "0")
data_2$Staging <- relevel(factor(data_2$Staging), ref = "3")  
data_2$Liver_metastasis <- relevel(factor(data_2$Liver_metastasis), ref = "0")
data_2$Lung_metastasis <- relevel(factor(data_2$Lung_metastasis), ref = "0")
data_2$Bone_metastasis <- relevel(factor(data_2$Bone_metastasis), ref = "0")
data_2$Treatment <- relevel(factor(data_2$Treatment), ref = "1")  
data_2$Line_of_therapy <- relevel(factor(data_2$Line_of_therapy), ref = "1")  
data_2$LDH <- relevel(factor(data_2$LDH), ref = "2")
data_2$CEA <- relevel(factor(data_2$CEA), ref = "2")
data_2$CA199 <- relevel(factor(data_2$CA199), ref = "2")
data_2$NSE <- relevel(factor(data_2$NSE), ref = "2")
data_2$Hb <- relevel(factor(data_2$Hb), ref = "2")
data_2$WBC <- relevel(factor(data_2$WBC), ref = "2")
data_2$AEC <- relevel(factor(data_2$AEC), ref = "2")
data_2$ANC <- relevel(factor(data_2$ANC), ref = "2")
data_2$ALC <- relevel(factor(data_2$ALC), ref = "2")
data_2$PLT <- relevel(factor(data_2$PLT), ref = "2")
data_2$NLR <- relevel(factor(data_2$NLR), ref = "2")
data_2$dNLR <- relevel(factor(data_2$dNLR), ref = "2")
data_2$PLR <- relevel(factor(data_2$PLR), ref = "2")
data_2$PNI <- relevel(factor(data_2$PNI), ref = "2")
data_2$SII <- relevel(factor(data_2$SII), ref = "2")
data_2$LWR <- relevel(factor(data_2$LWR), ref = "2")

# Initial Cox regression model (coxph)
cox_model_initial <- coxph(Surv(OS, Death) ~ age + sex + smoke + alcohol + Staging + Liver_metastasis + Lung_metastasis + Bone_metastasis + Line_of_therapy + Treatment + LDH + CEA + CA199 + NSE + Hb + WBC + AEC + ANC + ALC + PLT + NLR + dNLR + PLR + PNI + SII + LWR, 
                           data = data_2)
# Stepwise backward regression analysis
step_model <- step(cox_model_initial, direction = "backward")
# Calculating the C-index of model
concordance(step_model)$concordance - 1.96 * sqrt(concordance(step_model)$var)
concordance(step_model)$concordance + 1.96 * sqrt(concordance(step_model)$var)
c_index <- concordance(step_model)$concordance
print(paste("C-index:", c_index)) 

# ###
summary_step_model <- summary(step_model)
summary_df <- as.data.frame(summary_step_model$coefficients)
summary_df$exp_conf_low <- summary_step_model$conf.int[, "lower .95"]
summary_df$exp_conf_high <- summary_step_model$conf.int[, "upper .95"]
colnames(summary_df) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "exp_conf_low", "exp_conf_high")
summary_df <- cbind(Variable = rownames(summary_df), summary_df)
# ##
write_xlsx(summary_df, path = "C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2/cox_model_results.xlsx")






##############nomogram
#Data distribution settings
dd <- datadist(data_2)
options(datadist = "dd")
# Extracting variables for stepwise regression models
model_formula <- formula(step_model)

# Modify the variable name and corresponding value
data_2$Staging <- factor(data_2$Staging, levels = c(3, 4), labels = c("III", "IV"))
data_2$Bone_metastasis <- factor(data_2$Bone_metastasis, levels = c(0, 1), labels = c("No", "Yes"))
data_2$Line_of_therapy <- factor(data_2$Line_of_therapy, levels = c(1, 2), labels = c("First-line", "≥Second-line"))
data_2$Treatment <- factor(data_2$Treatment, levels = c(1, 2), labels = c("Immunotherapy/Immunotherapy+Target therapy", "Immunotherapy+Chemotherapy"))
data_2$LDH <- factor(data_2$LDH, levels = c(1, 2), labels = c(expression(">240"), expression("≤240")))
data_2$CEA <- factor(data_2$CEA, levels = c(1, 2), labels = c(expression(">5"), expression("≤5")))
data_2$CA199 <- factor(data_2$CA199, levels = c(1, 2), labels = c(expression(">37"), expression("≤37")))
data_2$ALC <- factor(data_2$ALC, levels = c(1, 2), labels = c(expression(">1.7"), expression("≤1.7")))
data_2$PNI <- factor(data_2$PNI, levels = c(1, 2), labels = c(expression(">49.5"), expression("≤49.5")))
data_2$SII <- factor(data_2$SII, levels = c(1, 2), labels = c(expression(">589.41"), expression("≤589.41")))

# Constructing the nomogram
nomogram_model <- cph(model_formula, data = data_2, x = TRUE, y = TRUE, surv = TRUE)

# Creating Survival Functions and Nomograms
surv <- Survival(nomogram_model)           
surv1 <- function(x) surv(12, lp = x)
surv2 <- function(x) surv(24, lp = x)
surv3 <- function(x) surv(36, lp = x)

nom <- rms::nomogram(
  nomogram_model,
  fun = list(surv1, surv2, surv3),
  lp = F,
  funlabel = c("1-Year Survival", "2-Year Survival", "3-Year Survival"),
  maxscale = 100,
  
  fun.at = c("0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1")
)

# ###
cairo_pdf("C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2/nomogram.pdf", width = 10, height = 7)
plot(nom, data= data_2, cex.var = 0.9, cex.axis = 0.7, lwd = 4)
dev.off()



#############Calibration curves
# ##
dd <- datadist(data_2)
options(datadist = "dd")

model_formula <- formula(step_model)

nomogram_model <- cph(model_formula, data = data_2, x = TRUE, y = TRUE, surv = TRUE)
# 1-year
calibration_1_year <- calibrate(nomogram_model, cmethod = "KM", method = "boot", 
                                u = 12, m = 100, B = 2000)
# 2-year
calibration_2_year <- calibrate(nomogram_model, cmethod = "KM", method = "boot", 
                                u = 24, m = 100, B = 2000)
# 3-year
calibration_3_year <- calibrate(nomogram_model, cmethod = "KM", method = "boot", 
                                u = 36, m = 100, B = 2000)
# plot
plot_calibration_curve <- function() {
  plot(calibration_1_year, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nomogram-predicted OS (%)", ylab = "Observed OS (%)", 
       lwd = 1, col = "palegreen3", sub = FALSE)
  plot(calibration_2_year, xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "", ylab = "", lwd = 1, col = "dodgerblue4", sub = FALSE, add = TRUE)
  plot(calibration_3_year, xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "", ylab = "", lwd = 1, col = "firebrick", sub = FALSE, add = TRUE)
  
  legend('bottomright', c('1-year', '2-year', '3-year'),
         col = c("palegreen3", "dodgerblue4", "firebrick"), lwd = 3, bty = 'n')
}

# save
cairo_pdf("C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2/calibration_curve.pdf", width = 8, height = 6)
plot_calibration_curve()
dev.off()








# Calculation of risk score
data_2$Predict_Score <- predict(nomogram_model, type = "lp")
# 使用三分位数分三组
data_2$Score_Group_Three <- cut(data_2$Predict_Score, 
                                breaks = quantile(data_2$Predict_Score, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                                include.lowest = TRUE, 
                                labels = c("Low", "Moderate", "High"))
# Calculate the three-quartile cutoff value
quantiles_3 <- quantile(data_2$Predict_Score, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
print(quantiles_3)


# save
write_xlsx(data_2, path = "C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2/Risk_groups.xlsx")

######### survival curves

fit_three <- survfit(Surv(OS, Death) ~ Score_Group_Three, data = data_2)
surv_plot <- ggsurvplot(
  fit_three, 
  data = data_2,
  conf.int = FALSE, 
  pval = TRUE, 
  risk.table = TRUE, 
  risk.table.height = 0.25,     
  palette = c("orange", "paleturquoise3","darkgreen"), 
  ggtheme = theme_classic() + theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.grid.major.y = element_line(color = "gray")),
  risk.table.y.text.col = TRUE, 
  risk.table.y.text = TRUE, 
  legend.title = "Risk Group",
  legend.labs = c("Low", "Moderate", "High"),
  title = "Overall survival based on Risk Group",
  xlab = "Time (months)",
  ylab = "Overall Survival (%)", 
  break.time.by = 12,   
  surv.scale = "percent"  
)

surv_plot$plot <- surv_plot$plot + scale_y_continuous(labels = c(0, 25, 50, 75, 100))
# save
cairo_pdf("OS_Three_group.pdf", width = 8, height = 6)
print(surv_plot)
dev.off()


# the median OS
median_OS <- surv_median(fit_three)
print(median_OS)
# HR（ "Low" is the reference group）
data_2$Score_Group_Three <- relevel(factor(data_2$Score_Group_Three), ref = "Low")

cox_model_three <- coxph(Surv(OS, Death) ~ Score_Group_Three, data = data_2)

summary(cox_model_three)


# Calculate 1-, 2-, and 3-year survival probabilities for each risk group
summary_fit <- summary(fit_three, times = c(12, 24, 36))
print(summary_fit)
# ##
risk_group_surv <- data.frame(
  Time = summary_fit$time,
  Risk_Group = rep(c("Low", "Moderate", "High"), each = 3),  # 根据分组
  Survival_Probability = summary_fit$surv,
  Lower_CI = summary_fit$lower,
  Upper_CI = summary_fit$upper
)
# save
write_xlsx(risk_group_surv, path = "survival_probability_in_risk_group.xlsx")







#######################AUC
# ##
surv_1_year <- predict(step_model, newdata = data_2, type = "survival", times = 12)
surv_2_year <- predict(step_model, newdata = data_2, type = "survival", times = 24)
surv_3_year <- predict(step_model, newdata = data_2, type = "survival", times = 36)

# ##
marker_1_year <- surv_1_year
marker_2_year <- surv_2_year
marker_3_year <- surv_3_year

# ROC curves
roc_1_year <- timeROC(T = data_2$OS, delta = data_2$Death, marker = marker_1_year,
                      cause = 1, weighting = "marginal", times = 12, iid = TRUE)
roc_2_year <- timeROC(T = data_2$OS, delta = data_2$Death, marker = marker_2_year,
                      cause = 1, weighting = "marginal", times = 24, iid = TRUE)
roc_3_year <- timeROC(T = data_2$OS, delta = data_2$Death, marker = marker_3_year,
                      cause = 1, weighting = "marginal", times = 36, iid = TRUE)

# se and 95% CI
se_1_year <- roc_1_year$inference$vect_sd_1[2]
se_2_year <- roc_2_year$inference$vect_sd_1[2]
se_3_year <- roc_3_year$inference$vect_sd_1[2]

auc_1_year <- roc_1_year$AUC[2]
auc_2_year <- roc_2_year$AUC[2]
auc_3_year <- roc_3_year$AUC[2]
auc_1_year_ci_lower <- auc_1_year - 1.96 * se_1_year
auc_1_year_ci_upper <- auc_1_year + 1.96 * se_1_year
auc_2_year_ci_lower <- auc_2_year - 1.96 * se_2_year
auc_2_year_ci_upper <- auc_2_year + 1.96 * se_2_year
auc_3_year_ci_lower <- auc_3_year - 1.96 * se_3_year
auc_3_year_ci_upper <- auc_3_year + 1.96 * se_3_year

# PDF
pdf("survival_ROC_1_2_3_years.pdf", width = 8, height = 6)

# Plotting ROC curves
plot(roc_1_year, time = 12, col = "orange", lty = 1, lwd = 2, xlab = "1-Specificity", ylab = "Sensitivity")
plot(roc_2_year, time = 24, col = "paleturquoise3", lty = 1, add = TRUE, lwd = 2)
plot(roc_3_year, time = 36, col = "darkgreen", lty = 1, add = TRUE, lwd = 2)

# 
legend("bottomright", legend = c(
  paste0("1-Year AUC (95% CI) = ", round(auc_1_year, 3), " (", round(auc_1_year_ci_lower, 3), " - ", round(auc_1_year_ci_upper, 3), ")"),
  paste0("2-Year AUC (95% CI) = ", round(auc_2_year, 3), " (", round(auc_2_year_ci_lower, 3), " - ", round(auc_2_year_ci_upper, 3), ")"),
  paste0("3-Year AUC (95% CI) = ", round(auc_3_year, 3), " (", round(auc_3_year_ci_lower, 3), " - ", round(auc_3_year_ci_upper, 3), ")")
), col = c("orange", "paleturquoise3","darkgreen"), lty = 1, lwd = 2)

dev.off()





##########Bootstrap validates C-index
# Bootstrap algorithm
boot_cox <- function(data, indices) {
  d <- data[indices, ]  # 进行bootstrap抽样
  fit <- coxph(formula(step_model), data = d)
  return(concordance(fit)$concordance)
}

# run Bootstrap 
set.seed(123)  
results <- boot(data = data_2, statistic = boot_cox, R = 1000)

#Check that 1000 samples conform to a normal distribution
hist(results$t, breaks = 30, main = "Bootstrap C-index Distribution", xlab = "C-index", col = "lightblue")
shapiro_test <- shapiro.test(results$t)
print(shapiro_test)

#Calculation of uncorrected_cindex and confidence intervals (mean)
uncorrected_cindex <- mean(results$t)
se_uncorrected <- sd(results$t) / sqrt(1000)  
ci_uncorrected_lower <- uncorrected_cindex - 1.96 * se_uncorrected
ci_uncorrected_upper <- uncorrected_cindex + 1.96 * se_uncorrected

# Calculation of Bias-corrected C-index and confidence intervals (BCA method)
bias_corrected_cindex <- 2 * c_index - uncorrected_cindex
ci_bca <- boot.ci(results, type = "bca")

# save
boot_data <- data.frame(
  Bootstrap_Replicate = 1:1000,
  Bootstrap_C_Index = results$t
)

ci_data <- data.frame(
  Statistic = c("Original C-index", "Uncorrected C-index", "Uncorrected C-index Lower Bound", "Uncorrected C-index Upper Bound", "Bias-corrected C-index", "Lower Bound", "Upper Bound"),
  Value = c(c_index, uncorrected_cindex, ci_uncorrected_lower, ci_uncorrected_upper, bias_corrected_cindex, ci_bca$bca[4], ci_bca$bca[5])
)
write_xlsx(list(Bootstrap_Results = boot_data, CI_Results = ci_data), 
           path = "C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2/bootstrap_c_index.xlsx")




##########Bootstrap validates 1/2/3 year AUCs
# 定义Bootstrap函数，收集每次Bootstrap的ROC曲线和AUC值
boot_roc_auc <- function(data, indices) {
  d <- data[indices, ]
  formula_step_model <- formula(step_model)    #重新设置step_model的公式环境
  environment(formula_step_model) <- environment()
  fit <- coxph(formula_step_model, data = d)
  
  # ##
  surv_prob_1_year <- predict(fit, newdata = d, type = "survival", times = 12)
  surv_prob_2_year <- predict(fit, newdata = d, type = "survival", times = 24)
  surv_prob_3_year <- predict(fit, newdata = d, type = "survival", times = 36)
  # Calculate ROC curve and AUC
  roc_1_year <- timeROC(T = d$OS, delta = d$Death, marker = surv_prob_1_year,
                        cause = 1, weighting = "marginal", times = 12, iid = TRUE)
  roc_2_year <- timeROC(T = d$OS, delta = d$Death, marker = surv_prob_2_year,
                        cause = 1, weighting = "marginal", times = 24, iid = TRUE)
  roc_3_year <- timeROC(T = d$OS, delta = d$Death, marker = surv_prob_3_year,
                        cause = 1, weighting = "marginal", times = 36, iid = TRUE)
  
  # ##
  return(c(roc_1_year$AUC[2], roc_2_year$AUC[2], roc_3_year$AUC[2]))
}    

# Bootstrap
set.seed(123)
results_roc_auc <- boot(data = data_2, statistic = boot_roc_auc, R = 1000)

# ##
auc_values <- results_roc_auc$t

#Check AUC distribution
hist(auc_values[, 1], main = "1-Year AUC Distribution", xlab = "AUC", col = "lightblue", breaks = 20)
hist(auc_values[, 2], main = "2-Year AUC Distribution", xlab = "AUC", col = "lightblue", breaks = 20)
hist(auc_values[, 3], main = "3-Year AUC Distribution", xlab = "AUC", col = "lightblue", breaks = 20)
# Shapiro-Wilk
shapiro_test_1_year <- shapiro.test(auc_values[, 1])
shapiro_test_2_year <- shapiro.test(auc_values[, 2])
shapiro_test_3_year <- shapiro.test(auc_values[, 3])
# 
print(shapiro_test_1_year)
print(shapiro_test_2_year)
print(shapiro_test_3_year)

# Calculation of AUC values (median) and 95% CI
mean_auc_1_year <- median(auc_values[, 1], na.rm = TRUE)
mean_auc_2_year <- median(auc_values[, 2], na.rm = TRUE)
mean_auc_3_year <- median(auc_values[, 3], na.rm = TRUE)
ci_auc_1_year <- boot.ci(results_roc_auc, type = "perc", index = 1)
ci_auc_2_year <- boot.ci(results_roc_auc, type = "perc", index = 2)
ci_auc_3_year <- boot.ci(results_roc_auc, type = "perc", index = 3)
# save
auc_results <- data.frame(
  Year = c("1-Year", "2-Year", "3-Year"),
  Mean_AUC = c(mean_auc_1_year, mean_auc_2_year, mean_auc_3_year),
  Lower_CI = c(ci_auc_1_year$percent[4], ci_auc_2_year$percent[4], ci_auc_3_year$percent[4]),
  Upper_CI = c(ci_auc_1_year$percent[5], ci_auc_2_year$percent[5], ci_auc_3_year$percent[5])
)
write_xlsx(auc_results, "C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/430_2/bootstrap_auc.xlsx")
