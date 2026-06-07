# ==============================================================================
# Simulation: Honest Selection-Free (Split-Group) vs. Double Dipping
# ==============================================================================

# 1. SETUP
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("energy")) install.packages("energy")
if (!require("furrr")) install.packages("furrr")

# ==============================================================================
# Definitive Reproducible Simulation: Honest Split-Group vs. Double Dipping
# ==============================================================================

# 1. SETUP
# استانداردسازی مولد اعداد تصادفی برای محیط‌های موازی (بسیار مهم)
RNGkind(kind = "L'Ecuyer-CMRG")

library(tidyverse)
library(energy)
library(furrr)


# بذر اصلی مقاله
MAIN_SEED <- 20250815
set.seed(MAIN_SEED)

# 2. CORE FUNCTIONS
get_disco_stat <- function(x_vals, g_labels, alpha_val) {
  # محاسبه آماره بدون پرموتیشن (فاز انتخاب)
  res <- energy::disco(x_vals, factors = g_labels, index = alpha_val, R = 0)
  return(res$statistic)
}

run_competition <- function(n, effect_type, dist_type) {
  # توجه: در این لایه Seed نباید به صورت دستی تنظیم شود چون توسط future_map مدیریت می‌شود
  
  # الف. تولید داده
  k <- 2
  if (effect_type == "null") {
    means <- c(0, 0); sds <- c(1, 1)
  } else {
    # یک اختلاف معنی‌دار برای واریانس (به همراه اختلاف در میانگین برای قدرت بیشتر)
    means <- c(0, 0.4); sds <- c(1, 1.8) 
  }
  
  if(dist_type == "normal") {
    x <- c(rnorm(n, means[1], sds[1]), rnorm(n, means[2], sds[2]))
  } else {
    x <- c(rt(n, df=3)*sds[1], rt(n, df=3)*sds[2])
  }
  g <- factor(rep(1:k, each=n))
  
  # فضای جستجوی گسترده (۲۰ کاندیدا برای ایجاد سوگیری شدید در Double Dipping)
  alpha_candidates <- seq(0.1, 2.0, by = 0.1)
  
  # --- استراتژی ۱: آزمون استاندارد ثابت (Alpha=1.0) ---
  p_fixed <- energy::disco(x, factors = g, index = 1.0, R = 499)$p.value # تعداد پرموتیشن بالاتر برای ثبات
  
  # --- استراتژی ۲: تقلب (Double Dipping) روی تمام داده ---
  all_ps <- sapply(alpha_candidates, function(a) {
    energy::disco(x, factors = g, index = a, R = 499)$p.value
  })
  p_double <- min(all_ps)
  
  # --- استراتژی ۳: پیشنهادی (Split-Group Honest) ---
  # استفاده از تقسیم‌بندی Stratified (برای تضمین حضور هر دو گروه در هر بخش)
  n_total <- length(x)
  grp1_idx <- which(g == "1")
  grp2_idx <- which(g == "2")
  
  setA_idx <- c(sample(grp1_idx, size = floor(n/2)), 
                sample(grp2_idx, size = floor(n/2)))
  
  x_A <- x[setA_idx]; g_A <- g[setA_idx]
  x_B <- x[-setA_idx]; g_B <- g[-setA_idx]
  
  # فاز انتخاب در نیمه A
  stats_A <- sapply(alpha_candidates, function(a) get_disco_stat(x_A, g_A, a))
  best_alpha <- alpha_candidates[which.max(stats_A)]
  
  # فاز استنتاج در نیمه B
  p_split <- energy::disco(x_B, factors = g_B, index = best_alpha, R = 499)$p.value
  
  return(data.frame(
    Strategy_Fixed = p_fixed,
    Strategy_DoubleDipping = p_double,
    Strategy_SplitGroup = p_split
  ))
}

# 3. EXECUTION
# در این بخش تعداد هسته‌ها و بذر فیکس می‌شود
plan(multisession, workers = parallel::detectCores() - 1)

# هر دو روی N=150 تنظیم شدند برای "موازنه قدرت" و "ثبات عددی"
cat("Running Final Analysis (n=150)... This will take a while.\n")

res_null <- future_map_dfr(1:400, ~run_competition(n=150, "null", "normal"), 
                           .options = furrr_options(seed = MAIN_SEED))

res_power <- future_map_dfr(1:400, ~run_competition(n=150, "power", "normal"), 
                            .options = furrr_options(seed = MAIN_SEED))
# 4. REPORT
cat("\n--- RESULT (N=150) ---\n")
list(NULL_RESULT = summarise(res_null, across(everything(), list(rate = ~mean(. < 0.05)))),
     POWER_RESULT = summarise(res_power, across(everything(), list(power = ~mean(. < 0.05))))) %>% print()