# VOC-Targeted--Breathomics--biomarker-pipeline
# ==============================================================================
# Script Name: VOC-Targeted-Breathomics-biomarker-pipeline.R
# Description: End-to-end pipeline for exhaled breath VOC metabolomics analysis.
#              Characterizes stable metabolic profiles in long-term patients.
# Features: ComplexHeatmap, 3D PCA, Modular Limma-LASSO-ML screening, ROCs, Boxplots
# Global Seed: 123
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup and Configuration
# ------------------------------------------------------------------------------
base_dir <- "F:/VOC/analysis/PN_VOC"
setwd(base_dir)

output_dir <- file.path(base_dir, "Integrated_Pipeline_Outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_packages <- c(
  "tidyverse", "caret", "limma", "glmnet", "e1071", "randomForest", "logistf", 
  "pROC", "ComplexHeatmap", "plotly", "ggpubr", "patchwork", "writexl", "circlize"
)
new_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]
if (length(new_packages) > 0) install.packages(new_packages, dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(caret); library(limma); library(glmnet)
  library(e1071); library(randomForest); library(logistf); library(pROC)
  library(ComplexHeatmap); library(plotly); library(ggpubr); library(patchwork)
  library(writexl); library(circlize)
})

set.seed(123)

# ------------------------------------------------------------------------------
# 2. Global Data Loading and Safe Cleaning
# ------------------------------------------------------------------------------
raw_data_file <- "PN_DE_data.csv"
if (!file.exists(raw_data_file)) stop("Raw data file is missing!")

dat <- read.csv(raw_data_file, check.names = FALSE) %>%
  mutate(Sample = as.character(Sample))

clinical_cols <- c("Sample", "Group", "Sex", "Age", "BMI", "Smoke_status", 
                   "Drink_status", "Duration_of_exposure", "Complications", 
                   "FVC_pct_pred", "FEV1_pct_pred", "FEV1_FVC_pct")

voc_cols_original <- setdiff(names(dat), clinical_cols)
safe_voc_cols <- make.names(voc_cols_original, unique = TRUE)
name_map <- tibble(Original_name = voc_cols_original, Safe_name = safe_voc_cols)
colnames(dat)[colnames(dat) %in% voc_cols_original] <- safe_voc_cols

dat <- dat %>%
  mutate(
    Group = factor(Group, levels = c("HC", "DE", "PN")),
    Complications = factor(Complications, levels = c(0, 1), labels = c("No", "Yes")),
    Smoke_status = factor(Smoke_status, levels = c(1, 2, 3), labels = c("Never", "Current", "Former")),
    Duration_of_exposure = as.numeric(Duration_of_exposure),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI)
  )

# ------------------------------------------------------------------------------
# 3. Global Data Overview: Heatmap and 3D PCA
# ------------------------------------------------------------------------------
# [PCA Analysis]
pca_res <- prcomp(dat[, safe_voc_cols], center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x) %>% mutate(Group = dat$Group, Sample = dat$Sample)

fig_pca <- plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group,
                   colors = c("HC" = "#4DBBD5FF", "DE" = "#00A087FF", "PN" = "#E64B35FF"),
                   marker = list(size = 5, opacity = 0.8), text = ~Sample) %>%
  layout(title = "3D PCA of Breath VOC Profiles", scene = list(
    xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")))
htmlwidgets::saveWidget(as_widget(fig_pca), file.path(output_dir, "Global_3D_PCA.html"))

# [Complex Heatmap]
mat_scaled <- scale(as.matrix(dat[, safe_voc_cols]))
rownames(mat_scaled) <- dat$Sample
colnames(mat_scaled) <- voc_cols_original

col_ha <- HeatmapAnnotation(
  Group = dat$Group,
  col = list(Group = c("HC" = "#4DBBD5FF", "DE" = "#00A087FF", "PN" = "#E64B35FF")),
  show_annotation_name = TRUE
)

pdf(file.path(output_dir, "Global_VOC_Heatmap.pdf"), width = 12, height = 8)
Heatmap(t(mat_scaled), name = "Z-score", top_annotation = col_ha, 
        col = colorRamp2(c(-2, 0, 2), c("#3C5488FF", "white", "#E64B35FF")),
        show_row_names = TRUE, show_column_names = FALSE,
        cluster_columns = TRUE, cluster_rows = TRUE,
        row_names_gp = gpar(fontsize = 8))
dev.off()

# ------------------------------------------------------------------------------
# 4. Automated Pipeline: Feature Screening, Modeling, and Visualization Wrapper
# ------------------------------------------------------------------------------
run_discovery_pipeline <- function(master_data, case_group, control_groups, out_prefix) {
  
  cat(sprintf("\n>>> Starting Module: %s vs %s <<<\n", case_group, paste(control_groups, collapse = "+")))
  
  sub_dir <- file.path(output_dir, out_prefix)
  dir.create(sub_dir, showWarnings = FALSE)
  
  # Prepare data structure
  sub_data <- master_data %>%
    filter(Group %in% c(case_group, control_groups)) %>%
    mutate(
      Outcome = if_else(Group == case_group, 1L, 0L),
      Target_Group = if_else(Outcome == 1, "Case", "Control"),
      Target_Group = factor(Target_Group, levels = c("Control", "Case"))
    )
  
  # Split dataset (70/30)
  train_idx <- createDataPartition(sub_data$Outcome, p = 0.70, list = FALSE)
  train_set <- sub_data[train_idx, ]
  valid_set <- sub_data[-train_idx, ]
  
  # Covariates and Preprocessing
  pre_voc <- preProcess(train_set[, safe_voc_cols], method = c("YeoJohnson", "center", "scale"))
  tr_voc_z <- predict(pre_voc, train_set[, safe_voc_cols])
  va_voc_z <- predict(pre_voc, valid_set[, safe_voc_cols])
  
  pre_cov <- preProcess(train_set[, c("Age", "BMI", "Duration_of_exposure")], method = c("center", "scale"))
  tr_cov_z <- predict(pre_cov, train_set[, c("Age", "BMI", "Duration_of_exposure")]) %>% rename(Age_z=Age, BMI_z=BMI, Dur_z=Duration_of_exposure)
  va_cov_z <- predict(pre_cov, valid_set[, c("Age", "BMI", "Duration_of_exposure")]) %>% rename(Age_z=Age, BMI_z=BMI, Dur_z=Duration_of_exposure)
  
  make_cov_matrix <- function(df, z_df) {
    tibble(
      Age_z = z_df$Age_z, BMI_z = z_df$BMI_z, Dur_z = z_df$Dur_z,
      Smk_cur = as.integer(df$Smoke_status == "Current"),
      Smk_for = as.integer(df$Smoke_status == "Former"),
      Comp = as.integer(df$Complications == "Yes")
    )
  }
  tr_cov <- make_cov_matrix(train_set, tr_cov_z)
  va_cov <- make_cov_matrix(valid_set, va_cov_z)
  
  # [1. Limma Screening]
  limma_dat <- bind_cols(train_set %>% select(Outcome), tr_cov)
  design <- model.matrix(~ 0 + as.factor(Outcome) + Age_z + BMI_z + Dur_z + Smk_cur + Smk_for + Comp, data = limma_dat)
  colnames(design)[1:2] <- c("Control", "Case")
  
  fit <- lmFit(t(as.matrix(tr_voc_z)), design)
  contrast <- makeContrasts(Case - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  
  limma_res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("Safe_name") %>%
    filter(adj.P.Val < 0.05) %>% pull(Safe_name)
  
  if(length(limma_res) == 0) limma_res <- topTable(fit2, coef = 1, number = 10, adjust.method = "BH") %>% rownames_to_column("Safe_name") %>% pull(Safe_name)
  
  # [2. LASSO Purification]
  x_lasso <- cbind(as.matrix(tr_voc_z[, limma_res]), as.matrix(tr_cov))
  cv_lasso <- cv.glmnet(x_lasso, train_set$Outcome, family = "binomial", alpha = 1, nfolds = 5)
  coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  final_vocs_safe <- rownames(coef_lasso)[rownames(coef_lasso) %in% limma_res & abs(coef_lasso[, 1]) > 0]
  
  if(length(final_vocs_safe) == 0) final_vocs_safe <- limma_res[1:5]
  
  # [3. ML Modeling]
  tr_x <- cbind(as.matrix(tr_voc_z[, final_vocs_safe, drop=FALSE]), as.matrix(tr_cov))
  va_x <- cbind(as.matrix(va_voc_z[, final_vocs_safe, drop=FALSE]), as.matrix(va_cov))
  
  m_svm <- svm(x = tr_x, y = train_set$Target_Group, kernel = "radial", probability = TRUE)
  p_svm_tr <- attr(predict(m_svm, tr_x, probability = TRUE), "probabilities")[, "Case"]
  p_svm_va <- attr(predict(m_svm, va_x, probability = TRUE), "probabilities")[, "Case"]
  
  m_rf <- randomForest(x = tr_x, y = train_set$Target_Group, ntree = 500)
  p_rf_tr <- predict(m_rf, tr_x, type = "prob")[, "Case"]
  p_rf_va <- predict(m_rf, va_x, type = "prob")[, "Case"]
  
  # [4. Plot and Export Results]
  pdf(file.path(sub_dir, paste0(out_prefix, "_ML_ROC.pdf")), width = 11, height = 5)
  par(mfrow = c(1, 2))
  plot.roc(train_set$Outcome, p_svm_tr, col="#E64B35FF", main="Training ROC", legacy.axes=TRUE)
  plot.roc(train_set$Outcome, p_rf_tr, col="#4DBBD5FF", add=TRUE)
  legend("bottomright", legend=c("SVM", "Random Forest"), col=c("#E64B35FF", "#4DBBD5FF"), lwd=2)
  
  plot.roc(valid_set$Outcome, p_svm_va, col="#E64B35FF", main="Validation ROC", legacy.axes=TRUE)
  plot.roc(valid_set$Outcome, p_rf_va, col="#4DBBD5FF", add=TRUE)
  dev.off()
  
  # Boxplot
  plot_df <- bind_rows(train_set %>% mutate(Split="Train"), valid_set %>% mutate(Split="Valid"))
  box_list <- lapply(final_vocs_safe, function(v) {
    ggplot(plot_df, aes(x = Group, y = .data[[v]], fill = Group)) +
      geom_boxplot(alpha=0.6, outlier.shape=NA) + geom_jitter(width=0.2, size=1, alpha=0.5) +
      facet_wrap(~Split) + theme_classic() + labs(title = name_map$Original_name[name_map$Safe_name == v], y = "Intensity") +
      scale_fill_manual(values = c("HC" = "#4DBBD5FF", "DE" = "#00A087FF", "PN" = "#E64B35FF")) +
      theme(legend.position="none", plot.title=element_text(hjust=0.5, face="bold"))
  })
  pdf(file.path(sub_dir, paste0(out_prefix, "_Feature_Boxplots.pdf")), width = 10, height = 8)
  print(wrap_plots(box_list))
  dev.off()
  
  cat(sprintf(">>> %s module completed. Number of features: %d <<<\n", out_prefix, length(final_vocs_safe)))
}

# ------------------------------------------------------------------------------
# 5. Execute Multi-Strategy Pipeline
# ------------------------------------------------------------------------------
run_discovery_pipeline(dat, "PN", c("HC", "DE"), "PN_vs_HC_DE")
run_discovery_pipeline(dat, "DE", "HC", "DE_vs_HC")
run_discovery_pipeline(dat, "PN", "HC", "PN_vs_HC")
run_discovery_pipeline(dat, "PN", "DE", "PN_vs_DE")

cat("\nCongratulations, the global personal analysis pipeline has finished executing! Results have been saved to the Integrated_Pipeline_Outputs directory.\n")
