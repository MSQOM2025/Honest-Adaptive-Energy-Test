# = an=============================================================================
#
#    R Script for Simulation Study of the Split-Group Energy Test (F_E,split)
#    DEFINITIVE FINAL VERSION
#
# ==============================================================================
# For absolute reproducibility of the entire script run
set.seed(20250815) # Using today's date (YYYYMMDD format)
#This is good practice, but not strictly necessary because of the `furrr_options(seed = TRUE)

# 1. SETUP
required_packages <- c("tidyverse", "future", "furrr", "energy")
not_installed <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(not_installed) > 0) { install.packages(not_installed) }

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(energy))
select <- dplyr::select


# 2. CONTROL PANEL
TEST_MODE <- TRUE # SET TO FALSE FOR PUBLICATION RUN
if (TEST_MODE) { SIM_REPLICATIONS <- 100; PERMUTATIONS <- 99 } else { SIM_REPLICATIONS <- 2000; PERMUTATIONS <- 999 }

# 3. CORE FUNCTIONS
generate_data <- function(k, n_per_group, means, sds, dist_type = "normal") {
  observations <- vector("list", k)
  for (i in 1:k) {
    if (dist_type == "normal") {
      observations[[i]] <- rnorm(n = n_per_group[i], mean = means[i], sd = sds[i])
    } else if (dist_type == "t") {
      observations[[i]] <- rt(n = n_per_group[i], df = 5) * sds[i] + means[i]
    }
  }
  tibble(observation = unlist(observations), group = factor(rep(1:k, times = n_per_group)))
}

# The novel hold-out permutation test
# --- REPLACE your old fe_split_test function with this new version ---

fe_split_test <- function(data, split_ratio) {
  # Step 1: Data Partitioning with a flexible ratio
  data_split <- data %>%
    group_by(group) %>%
    # Use the split_ratio to determine the size of Set A
    mutate(set = if_else(row_number() <= floor(n() * split_ratio), "A", "B")) %>%
    ungroup()
  
  data_A <- filter(data_split, set == "A")
  data_B <- filter(data_split, set == "B")
  
  # Check if splits are valid (at least k samples in each set)
  if (nrow(data_A) < length(unique(data$group)) || nrow(data_B) < length(unique(data$group))) return(NA_real_)
  
  # Step 2: Calculate Observed Statistic from Set A
  t_obs_result <- try(disco(data_A$observation, factors = data_A$group, R = 0), silent = TRUE)
  if (inherits(t_obs_result, "try-error")) return(NA_real_)
  t_obs <- t_obs_result$statistic
  
  # Step 3: Generate Null Distribution from Set B
  original_labels_B <- data_B$group
  t_perm <- numeric(PERMUTATIONS)
  for (i in 1:PERMUTATIONS) {
    perm_labels <- sample(original_labels_B)
    perm_result <- try(disco(data_B$observation, factors = perm_labels, R = 0), silent = TRUE)
    t_perm[i] <- if (inherits(perm_result, "try-error")) NA else perm_result$statistic
  }
  
  # Step 4: Calculate p-value
  p_value <- (sum(t_perm >= t_obs, na.rm = TRUE) + 1) / (PERMUTATIONS + 1)
  return(p_value)
}

# Wrapper for the simulation scenarios
# --- REPLACE your old run_scenario_fe function with this new version ---

run_scenario_fe <- function(params) {
  # This part is the same
  means <- rep(0, 4); sds <- rep(1, 4)
  if (params$effect_type == "mean_shift") { means <- c(0, 0, 0.5, 0.8) }
  if (params$effect_type == "variance_shift") { sds <- c(1, 2, 3, 4) }
  
  # Parallelize the replications
  p_values_list <- replicate(SIM_REPLICATIONS, {
    sim_data <- generate_data(k = 4, n_per_group = rep(params$n, 4), means = means, sds = sds, dist_type = params$dist_type)
    
    # --- The only change is here: Pass the split_ratio to the test ---
    fe_split_p_value <- fe_split_test(sim_data, split_ratio = params$split_ratio)
    
    # We still run the standard test if it's not a split-ratio experiment
    # For the split-ratio experiment, we only need the power of our method.
    if (params$split_ratio == 0.5) { # Only run the full test for the main comparison
        disco_full_p_value <- disco(sim_data$observation, factors = sim_data$group, R = PERMUTATIONS)$p.value
    } else {
        disco_full_p_value <- NA
    }
    
    c(fe_split = fe_split_p_value, disco_full = disco_full_p_value)

  }, simplify = "matrix")
  
  power_results <- tibble(test = rownames(p_values_list), value = rowMeans(p_values_list < 0.05, na.rm = TRUE))
  
  as_tibble(params) %>%
    bind_cols(power_results) %>%
    mutate(
      result_metric = if_else(params$effect_type == "null", "type1_error", "power"),
      n_sims = SIM_REPLICATIONS
    )
}

# --- REPLACE your old SETUP AND RUN STUDY block (Section 4) with this ---

# 4. SETUP AND RUN STUDY
num_cores <- 3
plan(multisession, workers = num_cores)
cat(paste("Parallel processing enabled on", num_cores, "cores.\n"))

# --- Define the original simulation grid ---
simulation_grid_main <- expand_grid(
  n = c(30, 60, 100),
  effect_type = c("null", "mean_shift", "variance_shift"),
  dist_type = c("normal", "t"),
  split_ratio = 0.5 # The main comparison uses a 50/50 split
)

# --- DEFINE THE NEW SPLIT-RATIO EXPERIMENT GRID ---
simulation_grid_split_ratio <- expand_grid(
  n = 100, # Use a fixed, representative sample size
  effect_type = "mean_shift", # Test with a clear signal
  dist_type = "normal",
  split_ratio = c(0.1, 0.3, 0.5, 0.7, 0.9) # The different ratios to test
)

# Combine both grids
full_simulation_grid <- bind_rows(simulation_grid_main, simulation_grid_split_ratio)


cat("Starting full simulation study...\n")
start_time <- Sys.time()
sim_results <- future_map(
  .x = 1:nrow(full_simulation_grid),
  ~ run_scenario_fe(params = as.list(full_simulation_grid[.x, ])),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE, packages = c("tidyverse", "energy"))
)
end_time <- Sys.time()
cat("Simulation study finished. Total time:", format(end_time - start_time), "\n")


# --- The SAVE RESULTS block remains the same (Section 5) ---
final_results_table <- bind_rows(sim_results)
print(final_results_table)
write_csv(final_results_table, "simulation_results_fe_split_final_EXTENDED.csv") # New file name

# --- The PLOTTING block (Figures 1-3) can remain the same. ---
# It will automatically work on the subset of the data with split_ratio = 0.5

# ==========================================================
#  NEW PLOT: Power vs. Split Ratio (The key new figure)
# ==========================================================
cat("\nGenerating new plot: Power vs. Split Ratio...\n")

split_ratio_plot_data <- final_results_table %>%
  filter(
    n == 100,
    effect_type == "mean_shift",
    dist_type == "normal",
    test == "fe_split"
  )

plot_split_ratio <- ggplot(split_ratio_plot_data, aes(x = split_ratio, y = value)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 4, color = "darkblue") +
  labs(
    title = "Figure X: Statistical Power as a Function of the Split Ratio",
    subtitle = "For F_E,split test with n=100 per group, k=4 (Mean Shift)",
    x = "Split Ratio (Proportion of data in Set A)",
    y = "Empirical Power"
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 14)

print(plot_split_ratio)
ggsave("plot_fe_split_power_vs_ratio.png", plot = plot_split_ratio, width = 9, height = 6, dpi = 300)
cat("Power vs. Split Ratio plot saved to plot_fe_split_power_vs_ratio.png\n")

cat("\nAnalysis complete. All results and new plots have been saved.\n")