# Given data are:
non_obese_non_exposed <- 0.03
non_obese_exposed <- non_obese_non_exposed * 1.2
obese_non_exposed <- non_obese_non_exposed * (2 / 3)
obese_exposed <- obese_non_exposed * 1.2

# Create a data frame to display incidence rates
incidence_table <- data.frame(
  Group = c("Non-Obese, Non-Exposed", "Non-Obese, Exposed", "Obese, Non-Exposed", "Obese, Exposed"),
  Incidence_Rate = c(non_obese_non_exposed, non_obese_exposed, obese_non_exposed, obese_exposed)
)
print(incidence_table)

# Calculating Odds and Odds Ratio for Non-Obese group as an example
odds_non_exposed <- non_obese_non_exposed / (1 - non_obese_non_exposed)
odds_exposed <- non_obese_exposed / (1 - non_obese_exposed)
odds_ratio <- odds_exposed / odds_non_exposed

# Displaying Odds Ratio and Relative Risk
cat("Odds Ratio (Non-Obese):", odds_ratio, "\n")
cat("Relative Risk (given): 1.2\n")
