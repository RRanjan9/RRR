# (a) Relative risk for city A (relaxed gun laws)
shootings_A <- 50
population_A <- 100000
risk_A <- shootings_A / population_A

shootings_B <- 10
population_B <- 100000
risk_B <- shootings_B / population_B

# Relative risk of gun violence in city A
relative_risk_A <- risk_A / risk_B
cat("Relative Risk of gun violence in city A:", relative_risk_A, "\n")

# (b) Relative risk for city B (strict gun laws)
relative_risk_B <- risk_B / risk_A
cat("Relative Risk of gun violence in city B:", relative_risk_B, "\n")
cat("This means the risk of gun violence in city B is 20% of the risk in city A, or 5 times lower.:\n")

# (c) Additional questions to consider:
cat("Questions to consider before concluding association:\n")
cat("1. Are the populations comparable in demographics and economic factors?\n")
cat("2. Could other factors (like law enforcement or economic conditions) influence gun violence rates?\n")
cat("3. Are the data collection periods the same for both cities?\n")
cat("4. Is the gun law the only difference between the cities?\n")
cat("5. How is gun violence measured (e.g., homicides, accidents, suicides)?\n")





# Data
cases_no_calcium <- 75
cases_use_calcium <- 25
non_cases_no_calcium <- 25
non_cases_use_calcium <- 75

# (b) Calculate odds of exposure in cases and non-cases
odds_cases <- cases_no_calcium / cases_use_calcium
odds_non_cases <- non_cases_no_calcium / non_cases_use_calcium

cat("Odds of exposure in cases:", odds_cases, "\n")
cat("Odds of exposure in non-cases:", odds_non_cases, "\n")

# (c) Calculate the odds ratio using cross-product
odds_ratio <- (cases_no_calcium * non_cases_use_calcium) / (cases_use_calcium * non_cases_no_calcium)
cat("Odds Ratio:", odds_ratio, "\n")

# (d) Compare prevalence ratio and odds ratio
prevalence_ratio <- (cases_no_calcium / 100) / (non_cases_no_calcium / 100)
cat("Prevalence Ratio:", prevalence_ratio, "\n")
cat("The odds ratio magnifies the difference compared to the prevalence ratio.\n")

#Key Difference:
cat("The prevalence ratio (3) compares relative exposure percentages, 
while the odds ratio (9) compares the likelihood of exposure between groups. 
The odds ratio magnifies differences, especially for common events, making it a stronger measure of association.")
