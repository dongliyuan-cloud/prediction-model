library(survival)
library(dplyr)
library(writexl)
library(survminer)
library(timeROC)

# Setting up a working directory and reading data
setwd("C:/Users/lenovo/Desktop/写论文/免疫治疗队列/预测模型/验证")   # your working_ directory
validation_data <- read.csv("184.csv", fileEncoding = "GB2312")


########### fit model
# Set reference levels (consistent with when building the model)
validation_data$Staging <- relevel(factor(validation_data$Staging), ref = "3")  
validation_data$Bone_metastasis <- relevel(factor(validation_data$Bone_metastasis), ref = "0")
validation_data$Line_of_therapy <- relevel(factor(validation_data$Line_of_therapy), ref = "1")  
validation_data$Treatment <- relevel(factor(validation_data$Treatment), ref = "2")  
validation_data$CEA <- relevel(factor(validation_data$CEA), ref = "2")
validation_data$CA199 <- relevel(factor(validation_data$CA199), ref = "2")
validation_data$LDH <- relevel(factor(validation_data$LDH), ref = "2")
validation_data$ALC <- relevel(factor(validation_data$ALC), ref = "2")
validation_data$PNI <- relevel(factor(validation_data$PNI), ref = "2")
validation_data$SII <- relevel(factor(validation_data$SII), ref = "2")

# fit model
cox_validation <- coxph(Surv(OS, Death) ~ Staging + Bone_metastasis + Line_of_therapy + 
                          Treatment + CEA + CA199 + LDH + ALC + PNI + SII, 
                        data = validation_data)
# Calculate C-index and confidence intervals
concordance_validation <- concordance(cox_validation)
c_index_validation <- concordance_validation$concordance
ci_lower_validation <- c_index_validation - 1.96 * sqrt(concordance_validation$var)
ci_upper_validation <- c_index_validation + 1.96 * sqrt(concordance_validation$var)
# ##
print(paste("C-index:", round(c_index_validation, 3), 
            "95% CI:", round(ci_lower_validation, 3), "-", round(ci_upper_validation, 3)))






############Predictive scoring, risk stratification
validation_data$Predict_Score <- predict(nomogram_model, newdata = validation_data, type = "lp")

# Grouping using tertiles of training cohort risk scores as cutoff values
validation_data$Score_Group_Three_val <- cut(validation_data$Predict_Score, 
                                         breaks = quantiles_3,  # 训练集的三分位数
                                         include.lowest = TRUE, 
                                         labels = c("Low", "Moderate", "High"))
# save
write_xlsx(validation_data, path = "validation_results_184.xlsx")

# survival curves
fit_three_val <- survfit(Surv(OS, Death) ~ Score_Group_Three_val, data = validation_data)
surv_plot_val <- ggsurvplot(
  fit_three_val, 
  data = validation_data,
  conf.int = FALSE, 
  pval = TRUE, 
  risk.table = TRUE, 
  risk.table.height = 0.25,     
  palette = c("orange", "paleturquoise3", "darkgreen"), 
  ggtheme = theme_classic() + theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.grid.major.y = element_line(color = "gray")),
  risk.table.y.text.col = TRUE,
  risk.table.y.text = TRUE,
  legend.title = "Risk Group",
  legend.labs = c("Low", "Moderate", "High"),
  title = "Overall survival based on Risk Group (Validation Set)",
  xlab = "Time (months)",
  ylab = "Overall Survival (%)", 
  break.time.by = 12,
  surv.scale = "percent"
)
# 
surv_plot_val$plot <- surv_plot_val$plot + scale_y_continuous(labels = c(0, 25, 50, 75, 100))
# save
cairo_pdf("validation_OS_Three_group_184.pdf", width = 8, height = 6)
print(surv_plot_val)
dev.off()

# median OS
median_OS <- surv_median(fit_three_val)
print(median_OS)
# HR
validation_data$Score_Group_Three_val <- relevel(factor(validation_data$Score_Group_Three_val), ref = "Low")
# 
cox_model_three_val <- coxph(Surv(OS, Death) ~ Score_Group_Three_val, data = validation_data)
# 
summary(cox_model_three_val)


# Calculate 1-, 2-, and 3-year survival probabilities for each risk group
summary_fit_val <- summary(fit_three_val, times = c(12, 24, 36))
print(summary_fit_val)





#########AUC
# 
surv_prob_1_year_val <- predict(cox_validation, newdata = validation_data, type = "survival", times = 12)
surv_prob_2_year_val <- predict(cox_validation, newdata = validation_data, type = "survival", times = 24)
surv_prob_3_year_val <- predict(cox_validation, newdata = validation_data, type = "survival", times = 36)

# 
marker_1_year_val <- surv_prob_1_year_val
marker_2_year_val <- surv_prob_2_year_val
marker_3_year_val <- surv_prob_3_year_val

# 
roc_1_year_val <- timeROC(T = validation_data$OS, delta = validation_data$Death, marker = marker_1_year_val,
                          cause = 1, weighting = "marginal", times = 12, iid = TRUE)
roc_2_year_val <- timeROC(T = validation_data$OS, delta = validation_data$Death, marker = marker_2_year_val,
                          cause = 1, weighting = "marginal", times = 24, iid = TRUE)
roc_3_year_val <- timeROC(T = validation_data$OS, delta = validation_data$Death, marker = marker_3_year_val,
                          cause = 1, weighting = "marginal", times = 36, iid = TRUE)

# 
pdf("validation_ROC_1_2_3_years_184.pdf", width = 8, height = 6)
plot(roc_1_year_val, time = 12, col = "orange", lty = 1, lwd = 2, xlab = "1-Specificity", ylab = "Sensitivity")
plot(roc_2_year_val, time = 24, col = "paleturquoise3", lty = 1, add = TRUE, lwd = 2)
plot(roc_3_year_val, time = 36, col = "darkgreen", lty = 1, add = TRUE, lwd = 2)

# 
se_1_year_val <- roc_1_year_val$inference$vect_sd_1[2]
se_2_year_val <- roc_2_year_val$inference$vect_sd_1[2]
se_3_year_val <- roc_3_year_val$inference$vect_sd_1[2]

auc_1_year_val <- roc_1_year_val$AUC[2]
auc_2_year_val <- roc_2_year_val$AUC[2]
auc_3_year_val <- roc_3_year_val$AUC[2]

ci_auc_1_year_val <- c(auc_1_year_val - 1.96 * se_1_year_val, auc_1_year_val + 1.96 * se_1_year_val)
ci_auc_2_year_val <- c(auc_2_year_val - 1.96 * se_2_year_val, auc_2_year_val + 1.96 * se_2_year_val)
ci_auc_3_year_val <- c(auc_3_year_val - 1.96 * se_3_year_val, auc_3_year_val + 1.96 * se_3_year_val)

# 
legend("bottomright", legend = c(
  paste0("1-Year AUC (95% CI) = ", round(auc_1_year_val, 3), " (", round(ci_auc_1_year_val[1], 3), " - ", round(ci_auc_1_year_val[2], 3), ")"),
  paste0("2-Year AUC (95% CI) = ", round(auc_2_year_val, 3), " (", round(ci_auc_2_year_val[1], 3), " - ", round(ci_auc_2_year_val[2], 3), ")"),
  paste0("3-Year AUC (95% CI) = ", round(auc_3_year_val, 3), " (", round(ci_auc_3_year_val[1], 3), " - ", round(ci_auc_3_year_val[2], 3), ")")
), col = c("orange", "paleturquoise3", "darkgreen"), lty = 1, lwd = 2)
dev.off()



# 
dd <- datadist(validation_data)
options(datadist = "dd")

# Prediction of validation cohort data using nomogram_model from training cohort
validation_data$Predict_Score_validate <- predict(nomogram_model, newdata = validation_data, type = "lp")

# 
cox_validation <- cph(Surv(OS, Death) ~ Predict_Score_validate, data = validation_data, x = TRUE, y = TRUE, surv = TRUE)

# calibration curves
calibration_1_year_val <- calibrate(cox_validation, cmethod = "KM", method = "boot", 
                                    u = 12, m = 50, B = 2000)
# 
calibration_2_year_val <- calibrate(cox_validation, cmethod = "KM", method = "boot", 
                                    u = 24, m = 50, B = 2000)
# 
calibration_3_year_val <- calibrate(cox_validation, cmethod = "KM", method = "boot", 
                                    u = 36, m = 50, B = 2000)
# 
plot_calibration_curve_validation <- function() {
  plot(calibration_1_year_val, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nomogram-predicted OS (%)", ylab = "Observed OS (%)", 
       lwd = 1, col = "palegreen3", sub = FALSE)
  plot(calibration_2_year_val, xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "", ylab = "", lwd = 1, col = "dodgerblue4", sub = FALSE, add = TRUE)
  plot(calibration_3_year_val, xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "", ylab = "", lwd = 1, col = "firebrick", sub = FALSE, add = TRUE)
  
  # 
  legend('bottomright', c('1-year', '2-year', '3-year'),
         col = c("palegreen3", "dodgerblue4", "firebrick"), lwd = 3, bty = 'n')
}

# save
cairo_pdf("calibration_curve_validation_184.pdf", width = 8, height = 6)
plot_calibration_curve_validation()
dev.off()


