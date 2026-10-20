# ==============================================================================
# Script: run_real_golub_multisplit.R
# Title: Multi-Split Sensitivity Analysis on High-Dimensional Golub Leukemia Data
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("multtest", update = FALSE, ask = FALSE)
suppressPackageStartupMessages({
  library(multtest)
  library(tidyverse)
  library(energy)
})

# 1. بارگذاری و آماده‌سازی داده‌های واقعی گلب
data(golub)
X <- t(golub) # 38 بیمار در سطرها، 3051 ژن در ستون‌ها
group <- factor(golub.cl, labels = c("ALL", "AML")) # ALL=27, AML=11

cat("Dataset loaded: N =", nrow(X), "patients, p =", ncol(X), "genes.\n")
print(table(group))

# 2. تنظیمات
ALPHA_GRID <- seq(0.1, 2.0, by = 0.1)
B_SPLITS   <- 50
R_PERM     <- 999
set.seed(110)

# تابع کمکی استخراج آماره DISCO بدون جایگشت
get_disco_stat <- function(mat, grp, a) {
  energy::disco(mat, factors = grp, index = a, R = 0)$statistic
}

# 3. اجرای آزمون استاندارد ثابت با Alpha = 1.0 روی کل داده (برای مقایسه)
p_fixed_full <- energy::disco(X, factors = group, index = 1.0, R = R_PERM)$p.value
cat("\nStandard DISCO (Fixed alpha=1.0, Full Data) P-value:", format.pval(p_fixed_full, eps=0.001), "\n")

# 4. اجرای پروتکل Multi-Split و بررسی حساسیت افراز (پاسخ به داور 1)
cat("\nRunning 50-Split Partitioning Protocol on Golub data...\n")

split_results <- vector("list", B_SPLITS)

for (b in 1:B_SPLITS) {
  # افراز تصادفی طبقه‌بندی‌شده
  idx_A <- integer(0)
  for (lvl in levels(group)) {
    grp_idx <- which(group == lvl)
    n_A <- floor(length(grp_idx) * 0.5)
    idx_A <- c(idx_A, sample(grp_idx, size = n_A))
  }
  
  X_A <- X[idx_A, , drop = FALSE];  grp_A <- factor(group[idx_A])
  X_B <- X[-idx_A, , drop = FALSE]; grp_B <- factor(group[-idx_A])
  
  # گام 1: انتخاب آلفا روی نمونه A
  stats_A <- sapply(ALPHA_GRID, function(a) get_disco_stat(X_A, grp_A, a))
  best_alpha <- ALPHA_GRID[which.max(stats_A)]
  
  # گام 2: آزمون جایگشتی کامل روی نمونه B با آلفای فیکس‌شده
  res_B <- energy::disco(X_B, factors = grp_B, index = best_alpha, R = R_PERM)
  
  split_results[[b]] <- tibble(
    split_id       = b,
    selected_alpha = best_alpha,
    p_val          = res_B$p.value
  )
}

df_splits <- bind_rows(split_results)

# 5. تجمیع کوشی پی‌وال‌ها (Cauchy Combination)
p_vec <- df_splits$p_val
p_clipped <- pmin(pmax(p_vec, 1e-5), 1 - 1e-5)
cauchy_stat <- mean(tan((0.5 - p_clipped) * pi))
p_multi <- as.numeric(0.5 - (atan(cauchy_stat) / pi))
p_multi <- max(0, min(1, p_multi))

# 6. گزارش تحلیلی برای درج در مقاله
cat("\n======================================================\n")
cat("          GOLUB MULTI-SPLIT SENSITIVITY REPORT        \n")
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