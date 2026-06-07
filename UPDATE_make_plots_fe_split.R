# ==========================================================
#
#         R Script to Generate All Plots for F_E,split Paper
#
# ==========================================================
set.seed(110)
# --- 1. SETUP: Load libraries and the definitive simulation results ---
library(tidyverse)

# Load the results from the long simulation run.
# Make sure this CSV file is in your working directory.
final_results_table <- read_csv("E:/Salehi/question nad answers about AA new f test paper/Idea 3 F_E test/Code for F_E_Split/simulation_results_fe_split_final.csv")

# --- 2. DATA PREPARATION ---
plot_data <- final_results_table %>%
  # Make the test names more descriptive
  mutate(test_label = factor(case_when(
    test == "fe_split"   ~ "F_E,split (Ours)",
    # This needs to match the column name in your CSV, which was `disco_full.p-value`
    test == "disco_full.p-value" ~ "Standard DISCO Test" 
  ), levels = c("Standard DISCO Test", "F_E,split (Ours)"))) %>%
  # Create clean labels
  mutate(dist_label = factor(dist_type, levels=c("normal", "t"), labels=c("Normal Dist.", "t-Dist. (df=5)")))

# Calculate the ratio of powers
power_ratio_data <- plot_data %>%
  filter(result_metric == "power") %>%
  
  # Select only the columns we need for this plot
  select(n, effect_type, dist_label, test_label, value) %>%
  
  # Pivot the data. Critically, we now specify the 'id_cols' that
  # define a unique scenario.
  pivot_wider(
    names_from = test_label, 
    values_from = value,
    id_cols = c(n, effect_type, dist_label) # This tells pivot how to match rows
  ) %>%
  
  # Now the mutate will work because both columns will be populated.
  mutate(power_ratio = `F_E,split (Ours)` / `Standard DISCO Test`) %>%
  
  # Make effect_type labels cleaner for the plot
  mutate(effect_label = case_when(
    effect_type == "mean_shift" ~ "Mean Shift",
    effect_type == "variance_shift" ~ "Variance Shift"
  ))

# The ggplot code below this remains the same and will now work.


plot4 <- ggplot(power_ratio_data, aes(x = factor(n), y = power_ratio, group = dist_label, color = dist_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black") +
  facet_wrap(~effect_label) +
  labs(
    title = "Figure 4: Relative Power of the F_E,split Test",
    subtitle = "Ratio of F_E,split Power to Standard DISCO Power. The 'Price of Purity' diminishes as n grows.",
    x = "Sample Size per Group (n)",
    y = "Relative Power (Ratio)",
    color = "Distribution"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0.5, 1.05)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

print(plot4)
ggsave("plot_fe_split_relative_power.png", plot = plot4, width = 10, height = 6, dpi = 300)
cat("Relative Power plot saved to plot_fe_split_relative_power.png\n")cat("All plots have been generated and saved successfully.\n")