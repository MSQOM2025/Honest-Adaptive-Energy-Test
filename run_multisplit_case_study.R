# ==============================================================================
# Script: run_multisplit_case_study.R
# Title: Multi-Split Partitioning & Sensitivity Analysis (k=3 Multi-Group Setting)
# Author: Dr. Mohsen Salehi (Revision 2026)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(energy)
})

# ۱. تنظیم دانه اصلی (دقیقاً ۱۱۰ طبق فرمایش شما)
MAIN_SEED <- 110
set.seed(MAIN_SEED)

ALPHA_GRID <- seq(0.1, 2.0, by = 0.1) # بیست کاندیدا
B_SPLITS   <- 50                      # ۵۰ بار تکرار افراز تصادفی
R_PERM     <- 999                     # تعداد جایگشت‌ها

# ۲. بارگذاری داده‌ها (اگر از قبل در محیط لود شده باشند استفاده می‌شوند؛ وگرنه ساخته می‌شوند)
if (!exists("data_matrix") || !exists("group_labels")) {
  cat("--- Generating Case Study Dataset (N=150, p=100, k=3) ---\n")
  set.seed(12345)
  n_features <- 100
  n_samples_per_group <- 50
  
  group_a <- matrix(rnorm(n_samples_per_group * n_features, mean = 1.0, sd = 0.5), nrow = n_samples_per_group)
  group_b <- matrix(rnorm(n_samples_per_group * n_features, mean = 1.2, sd = 1.0), nrow = n_samples_per_group)
  group_c <- matrix(rnorm(n_samples_per_group * n_features, mean = 1.1, sd = 1.5), nrow = n_samples_per_group)
  
  data_matrix  <- rbind(group_a, group_b, group_c)
  group_labels <- factor(rep(c("Control", "Treated", "Non-responder"), each = n_samples_per_group))
} else {
  cat("--- Using Loaded Dataset from Database/Environment ---\n")
  data_matrix  <- as.matrix(data_matrix)
  group_labels <- factor(group_labels)
}

cat("Dataset Dimensions: N =", nrow(data_matrix), "subjects, p =", ncol(data_matrix), "features.\n")
cat("Sample size per group:\n")
print(table(group_labels))

# تابع استخراج آماره DISCO بدون جایگشت (برای مرحله Discovery)
get_disco_stat <- function(mat, grp, a) {
  energy::disco(mat, factors = grp, index = a, R = 0)$statistic
}

# ۳. اجرای آزمون استاندارد روی کل داده با Alpha = 1.0 (برای مقایسه)
cat("\nRunning Standard DISCO (Fixed alpha=1.0) on full sample...\n")
p_fixed_full <- energy::disco(data_matrix, factors = group_labels, index = 1.0, R = R_PERM)$p.value
cat("Standard DISCO (Fixed alpha=1.0) P-value:", format.pval(p_fixed_full, eps=0.001), "\n")

# ۴. اجرای پروتکل Multi-Split و بررسی حساسیت افراز (۵۰ بار تقسیم تصادفی)
cat("\nRunning 50-Split Partitioning Protocol...\n")
set.seed(MAIN_SEED) # تضمین تکرارپذیری در افرازها

split_results <- vector("list", B_SPLITS)

for (b in 1:B_SPLITS) {
  # افراز تصادفی طبقه‌بندی‌شده ۵۰/۵۰
  idx_A <- integer(0)
  for (lvl in levels(group_labels)) {
    grp_idx <- which(group_labels == lvl)
    n_A <- floor(length(grp_idx) * 0.5)
    idx_A <- c(idx_A, sample(grp_idx, size = n_A))
  }
  
  mat_A <- data_matrix[idx_A, , drop = FALSE];  grp_A <- factor(group_labels[idx_A])
  mat_B <- data_matrix[-idx_A, , drop = FALSE]; grp_B <- factor(group_labels[-idx_A])
  
  # گام ۱: انتخاب آلفای بهینه روی نمونه A
  stats_A <- sapply(ALPHA_GRID, function(a) get_disco_stat(mat_A, grp_A, a))
  best_alpha <- ALPHA_GRID[which.max(stats_A)]
  
  # گام ۲: آزمون جایگشتی مستقل روی نمونه B با آلفای انتخاب‌شده
  res_B <- energy::disco(mat_B, factors = grp_B, index = best_alpha, R = R_PERM)
  
  split_results[[b]] <- tibble(
    split_id       = b,
    selected_alpha = best_alpha,
    p_val          = res_B$p.value
  )
}

df_splits <- bind_rows(split_results)

# ۵. تجمیع کوشی پی‌وال‌ها (Cauchy Combination)
p_vec <- df_splits$p_val
p_clipped <- pmin(pmax(p_vec, 1e-5), 1 - 1e-5)
cauchy_stat <- mean(tan((0.5 - p_clipped) * pi))
p_multi <- as.numeric(0.5 - (atan(cauchy_stat) / pi))
p_multi <- max(0, min(1, p_multi))

# ۶. گزارش تحلیلی برای درج در مقاله
cat("\n======================================================\n")
cat("       MULTI-GROUP MULTI-SPLIT SENSITIVITY REPORT     \n")
cat("======================================================\n")
cat("Selected Alpha Distribution across 50 splits:\n")
print(table(df_splits$selected_alpha))

cat("\nP-value Summary across 50 random splits:\n")
cat(sprintf("Min P-value:    %.4f\n", min(df_splits$p_val)))
cat(sprintf("25th Percentile:%.4f\n", quantile(df_splits$p_val, 0.25)))
cat(sprintf("Median P-value: %.4f\n", median(df_splits$p_val)))
cat(sprintf("75th Percentile:%.4f\n", quantile(df_splits$p_val, 0.75)))
cat(sprintf("Max P-value:    %.4f\n", max(df_splits$p_val)))
cat(sprintf("\n>>> Aggregated Cauchy P-value (p_multi): %.4f <<<\n", p_multi))
cat("======================================================\n")