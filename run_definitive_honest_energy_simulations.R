# ==============================================================================
# Script: run_definitive_simulations.R
# Title: Definitive Multivariate Simulation Framework for Honest Adaptive DISCO
# Author: Dr. Mohsen Salehi (Revision 2026)
# ==============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 300)
# 1. SETUP AND DEPENDENCIES
required_pkgs <- c("dplyr", "tidyr","Rcpp","gsl","readr",  "tidyverse", "energy", "mvtnorm", "furrr", "future", "parallel")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(tidyverse)
  library(energy)
  library(mvtnorm)
  library(furrr)
  library(future)
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(Rcpp)
  library(gsl)
library(readr)

})

#install.packages("dplyr", dependencies = TRUE)
# تنظیم موتور تولید عدد تصادفی برای پردازش موازی دقیق
RNGkind(kind = "L'Ecuyer-CMRG")
MAIN_SEED <- 110
set.seed(MAIN_SEED)

# ==============================================================================
# 2. CONTROL PANEL (تنظیمات حالت تست یا اجرای نهایی)
# ==============================================================================
TEST_MODE <- TRUE  # اگر TRUE باشد ظرف چند دقیقه اجرا می‌شود؛ برای اجرای نهایی مقاله FALSE کنید

if (TEST_MODE) {
  SIM_REPS     <- 50   # تعداد تکرار مونت‌کارلو برای تست
  N_PERM       <- 99   # تعداد جایگشت‌ها
  B_SPLITS     <- 15   # تعداد تقسیم‌ها در روش Multi-Split
  SAMPLE_SIZES <- c(30, 60)
  DIM_LIST     <- c(10)
} else {
  SIM_REPS     <- 1000 # تعداد تکرار مونت‌کارلو برای نسخه نهایی ژورنال
  N_PERM       <- 500  # تعداد استاندارد جایگشت‌ها
  B_SPLITS     <- 50   # تعداد تقسیم‌ها طبق مقاله مینهوزن
  SAMPLE_SIZES <- c(30, 60, 100)
  DIM_LIST     <- c(10, 50)
}

ALPHA_GRID <- seq(0.1, 2.0, by = 0.1) # شبکه کاندیداهای آلفا (M=20)

# ==============================================================================
# 3. MULTIVARIATE DATA GENERATOR (تولید داده‌های چندمتغیره در بعد d)
# ==============================================================================
generate_multivariate_data <- function(k = 4, n = 30, d = 10, 
                                       dgm = c("gaussian_null", "t3_null", "gaussian_shift", "cauchy_shift"), 
                                       delta = 0.35) {
  dgm <- match.arg(dgm)
  data_list <- vector("list", k)
  
  for (j in 1:k) {
    # تنظیم بردار مکان (میانگین) در بعد d
    if (dgm %in% c("gaussian_null", "t3_null")) {
      mu_vec <- rep(0, d)
    } else {
      # گروه 1 در مبدا، گروه‌های بعدی دارای انتقال به میزان delta
      mu_vec <- if (j == 1) rep(0, d) else rep(delta, d)
    }
    
    # تولید ماتریس n x d بر اساس توزیع مورد نظر
    if (dgm == "gaussian_null" || dgm == "gaussian_shift") {
      X_j <- rmvnorm(n = n, mean = mu_vec, sigma = diag(d))
    } else if (dgm == "t3_null") {
      X_j <- rmvt(n = n, sigma = diag(d), df = 3) + matrix(mu_vec, nrow = n, ncol = d, byrow = TRUE)
    } else if (dgm == "cauchy_shift") {
      # توزیع تی با 1 درجه آزادی همان کوشی چندمتغیره (دم بی‌نهایت سنگین) است
      X_j <- rmvt(n = n, sigma = diag(d), df = 1) + matrix(mu_vec, nrow = n, ncol = d, byrow = TRUE)
    }
    data_list[[j]] <- X_j
  }
  
  X <- do.call(rbind, data_list)
  group <- factor(rep(1:k, each = n))
  return(list(X = X, group = group))
}

# ==============================================================================
# 4. CORE TEST ENGINES
# ==============================================================================

# تابع کمکی محاسبه آماره DISCO بدون جایگشت (برای مرحله Discovery)
get_disco_stat <- function(X, group, alpha_val) {
  res <- energy::disco(X, factors = group, index = alpha_val, R = 0)
  return(res$statistic)
}

# پروتکل تک‌تقسیمی صادقانه (F_E,split)
run_fe_split_once <- function(X, group, alpha_grid, split_ratio = 0.5, R_perm = 499) {
  # افراز تصادفی طبقه‌بندی‌شده (Stratified Partitioning)
  idx_A <- integer(0)
  for (lvl in levels(group)) {
    grp_idx <- which(group == lvl)
    n_A <- floor(length(grp_idx) * split_ratio)
    idx_A <- c(idx_A, sample(grp_idx, size = n_A))
  }
  
  X_A <- X[idx_A, , drop = FALSE];     group_A <- factor(group[idx_A])
  X_B <- X[-idx_A, , drop = FALSE];    group_B <- factor(group[-idx_A])
  
  # گام 1: انتخاب آلفای بهینه منحصراً روی نمونه A
  stats_A <- sapply(alpha_grid, function(a) get_disco_stat(X_A, group_A, a))
  best_alpha <- alpha_grid[which.max(stats_A)]
  
  # گام 2: آزمون جایگشتی مستقل روی نمونه B با آلفای فیکس‌شده
  p_val <- energy::disco(X_B, factors = group_B, index = best_alpha, R = R_perm)$p.value
  return(list(p_val = p_val, alpha = best_alpha))
}

# ارزیابی یک سناریوی کامل با تمام ۴ استراتژی
evaluate_single_replicate <- function(k, n, d, dgm, delta, alpha_grid, R_perm, B_splits) {
  dat <- generate_multivariate_data(k = k, n = n, d = d, dgm = dgm, delta = delta)
  X <- dat$X
  group <- dat$group
  
  # -------------------------------------------------------------
  # استراتژی ۱: آزمون استاندارد ثابت (Alpha = 1.0 a priori)
  # -------------------------------------------------------------
  p_fixed <- energy::disco(X, factors = group, index = 1.0, R = R_perm)$p.value
  
  # -------------------------------------------------------------
  # استراتژی ۲: انتخاب ناایو روی کل نمونه (Naive Full-Sample)
  # -------------------------------------------------------------
  # یافتن آلفایی که آماره کل داده را ماکزیمم کند
  full_stats <- sapply(alpha_grid, function(a) get_disco_stat(X, group, a))
  best_alpha_full <- alpha_grid[which.max(full_stats)]
  # تست با همان آلفا روی کل داده (سوگیری شدید انتخاب)
  p_naive <- energy::disco(X, factors = group, index = best_alpha_full, R = R_perm)$p.value
  
  # -------------------------------------------------------------
  # استراتژی ۳: تک‌تقسیمی صادقانه (F_E,split)
  # -------------------------------------------------------------
  res_single <- run_fe_split_once(X, group, alpha_grid, split_ratio = 0.5, R_perm = R_perm)
  p_split <- res_single$p_val
  
  # -------------------------------------------------------------
  # استراتژی ۴: چندتقسیمی تجمیعی (F_E,multi با فرمول Meinshausen)
  # -------------------------------------------------------------
  p_vec <- replicate(B_splits, {
    run_fe_split_once(X, group, alpha_grid, split_ratio = 0.5, R_perm = R_perm)$p_val
  })
  
  # تجمیع میانه (Quantile Aggregation با تتا = 0.5)
  # تجمیع استاندارد کوشی (Liu & Xie, 2020 JASA) که خطای نوع اول را دقیقاً روی 0.05 نگه می‌دارد:
   p_clipped <- pmin(pmax(p_vec, 1e-5), 1 - 1e-5)
   cauchy_stat <- mean(tan((0.5 - p_clipped) * pi))
   p_multi <- as.numeric(0.5 - (atan(cauchy_stat) / pi))
   p_multi <- max(0, min(1, p_multi))  
  return(c(
    Fixed_Alpha1 = p_fixed,
    Naive_Full   = p_naive,
    FE_Split     = p_split,
    FE_Multi     = p_multi
  ))
}

# ==============================================================================
# 5. EXECUTION PIPELINE (پردازش موازی)
# ==============================================================================
run_full_benchmark <- function() {
  num_cores <- max(1, parallel::detectCores() - 1)
  plan(multisession, workers = num_cores)
  cat(sprintf(">>> Parallel engine enabled on %d cores.\n", num_cores))
  cat(sprintf(">>> TEST_MODE = %s (Reps=%d, Permutations=%d, Splits=%d)\n", 
              TEST_MODE, SIM_REPS, N_PERM, B_SPLITS))
  
  # سناریوهای مورد نظر
  scenarios <- expand_grid(
    d = DIM_LIST,
    n = SAMPLE_SIZES,
    dgm = c("gaussian_null", "t3_null", "gaussian_shift", "cauchy_shift")
  ) %>%
    mutate(delta = case_when(
      dgm %in% c("gaussian_null", "t3_null") ~ 0.0,
      dgm == "gaussian_shift" ~ 0.30,
      dgm == "cauchy_shift"   ~ 0.60
    ))
  
  all_results <- vector("list", nrow(scenarios))
  
  for (s in 1:nrow(scenarios)) {
    sc <- scenarios[s, ]
    cat(sprintf("\n[%d/%d] Running scenario: d=%d, n=%d, DGM=%s, delta=%.2f ...\n",
                s, nrow(scenarios), sc$d, sc$n, sc$dgm, sc$delta))
    
    start_t <- Sys.time()
    
    # اجرای موازی تکرارها
    sim_out <- future_map_dfr(1:SIM_REPS, function(i) {
      p_vals <- evaluate_single_replicate(
        k = 4, n = sc$n, d = sc$d, dgm = sc$dgm, delta = sc$delta,
        alpha_grid = ALPHA_GRID, R_perm = N_PERM, B_splits = B_SPLITS
      )
      as.data.frame(t(p_vals))
    }, .options = furrr_options(seed = TRUE, packages = c("energy", "mvtnorm")))
    
    # محاسبه نرخ رد در سطح 0.05
    rejection_rates <- sim_out %>%
      summarise(across(everything(), ~ mean(. <= 0.05, na.rm = TRUE))) %>%
      mutate(
        d = sc$d,
        n = sc$n,
        dgm = sc$dgm,
        delta = sc$delta,
        sim_reps = SIM_REPS,
        type = if_else(grepl("null", sc$dgm), "Type1_Error", "Power")
      )
    
    all_results[[s]] <- rejection_rates
    elapsed <- round(difftime(Sys.time(), start_t, units = "secs"), 1)
    cat(sprintf("Completed in %s seconds.\n", elapsed))
    print(rejection_rates %>% select(d, n, dgm, type, everything()))
  }
  
  final_df <- bind_rows(all_results)
  write_csv(final_df, "definitive_simulation_results.csv")
  cat("\n>>> All scenarios completed successfully! Results saved to definitive_simulation_results.csv\n")
  return(final_df)
}

# اجرای شبیه‌سازی
simulation_start <- Sys.time()

cat("\n====================================\n")
cat("Simulation started at:",
    format(simulation_start, "%Y-%m-%d %H:%M:%S"),
    "\n")
cat("====================================\n")

results_table <- run_full_benchmark()


simulation_end <- Sys.time()

total_time <- difftime(
  simulation_end,
  simulation_start,
  units = "secs"
)

cat("\n====================================\n")
cat("Simulation finished at:",
    format(simulation_end, "%Y-%m-%d %H:%M:%S"),
    "\n")
cat("Total execution time:",
    round(as.numeric(total_time), 2),
    "seconds\n")
cat("Total execution time:",
    round(as.numeric(total_time) / 60, 2),
    "minutes\n")
cat("====================================\n")