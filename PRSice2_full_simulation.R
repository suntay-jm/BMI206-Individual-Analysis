# seed for reproducibility
set.seed(42)

# setting number of people, total number of SNPs, and causal SNPs
n_individuals = 1000
n_snps = 10000
causal_snps = 500  

# simulating gwas summary stats
simulate_gwas_stats = function(n_snps, causal_snps, beta_mean, beta_sd, noise_sd) {
  # randomly choosing which SNPs are causal 
  causal_indices = sample(1:n_snps, causal_snps)
  
  # simulate effect sizes (beta)
  beta = rep(0, n_snps) # starting with no effect for all SNPs
  beta[causal_indices] = rnorm(causal_snps, beta_mean, beta_sd) # adding effects for causal SNPs
  
  # simulate standard errors (SE)
  se = abs(rnorm(n_snps, mean = 0, sd = noise_sd))
  
  # calculate z-scores to measure each SNP's effect
  z_scores = beta / se
  
  # using z-scores to calculate p-vals for signficance
  p_values = 2 * pnorm(-abs(z_scores))  # two-tailed p-value for significance
  
  # return summary stats with SNP names, effects, SE, and p-vals
  data.frame(SNP = paste0("rs", 1:n_snps), Beta = beta, SE = se, P = p_values)
}

# generate GWAS summary stats for each population using simulate_gwas_stats()
gwas_eur = simulate_gwas_stats(n_snps, causal_snps, beta_mean = 0.2, beta_sd = 0.05, noise_sd = 0.1)
gwas_afr = simulate_gwas_stats(n_snps, causal_snps, beta_mean = 0.05, beta_sd = 0.02, noise_sd = 0.1)
gwas_eas = simulate_gwas_stats(n_snps, causal_snps, beta_mean = 0.1, beta_sd = 0.03, noise_sd = 0.1)

# generating genotype matrices
generate_genotypes = function(n_individuals, n_snps, maf) {
  
  # giving each person 0, 1, or 2 copies of a certain SNP based on MAF
  matrix(rbinom(n_individuals * n_snps, 2, maf), nrow = n_individuals, ncol = n_snps)
}

# MAFs for each population
maf_eur = 0.2
maf_afr = 0.3
maf_eas = 0.1

genotypes_eur = generate_genotypes(n_individuals, n_snps, maf_eur)
genotypes_afr = generate_genotypes(n_individuals, n_snps, maf_afr)
genotypes_eas = generate_genotypes(n_individuals, n_snps, maf_eas)

# writing simulating phenotype function (y = X * Beta + e)
simulate_phenotypes = function(genotypes, gwas_stats, noise_sd) {
  beta = gwas_stats$Beta
  e = rnorm(nrow(genotypes), mean = 0, sd = noise_sd)
  y = genotypes %*% beta + e
  list(y = y, beta = beta)
}

# using simulate_phenotypes on each population
phenotypes_eur = simulate_phenotypes(genotypes_eur, gwas_eur, noise_sd = 0.5)
phenotypes_afr = simulate_phenotypes(genotypes_afr, gwas_afr, noise_sd = 0.5)
phenotypes_eas = simulate_phenotypes(genotypes_eas, gwas_eas, noise_sd = 0.5)

# function to create mock PRS
generate_mock_prs = function(population, max_r2, max_p_threshold) {
  thresholds = seq(0.01, 1, by = 0.01) # testing different p-val thresholds
  r2_values = dnorm(thresholds, mean = max_p_threshold, sd = 0.1) * max_r2 # R^2 depends on threshold
  r2_values[r2_values < 0] = 0  #  no negative R^2 values
  data.frame(Population = population, P_Threshold = thresholds, R2 = r2_values)
}

# simulated PRS results 
mock_eur = generate_mock_prs("EUR", max_r2 = 0.25, max_p_threshold = 0.05)
mock_afr = generate_mock_prs("AFR", max_r2 = 0.1, max_p_threshold = 0.2)
mock_eas = generate_mock_prs("EAS", max_r2 = 0.15, max_p_threshold = 0.1)

# visualizing
library(ggplot2)

plot_prs = function(mock_data, title_color) {
  ggplot(mock_data, aes(x = P_Threshold, y = R2)) +
    geom_bar(stat = "identity", width = 0.02, fill = title_color) +
    geom_bar(data = mock_data %>% filter(R2 == max(R2)),  # highlight max R^2
             stat = "identity", width = 0.02, fill = NA, color = "black", size = 1.2) +
    geom_text(data = mock_data %>% filter(R2 == max(R2)),  # labeling max R^2 
              aes(label = paste0("R² = ", round(R2, 4), ", P_T = ", round(P_Threshold, 4))),
              vjust = 0.5, hjust = -0.1, size = 4) +
    labs(title = paste("PRS Performance:", mock_data$Population[1]),
         x = "P-value Threshold (P_T)",
         y = expression(R^2)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
}

# plot EUR, AFR, and EAS with specific colors
plot_prs(mock_eur, "cornflowerblue")
plot_prs(mock_afr, "burlywood")
plot_prs(mock_eas, "plum3")



# sensitivity analysis
sensitivity_analysis = function(mock_data) {
  # find the row with the maximum R2
  max_r2_row = mock_data %>% filter(R2 == max(R2))
  max_threshold = max_r2_row$P_Threshold
  
  # create Threshold_1 as the fixed max R2 threshold
  fixed_threshold = rep(max_threshold, nrow(mock_data))
  
  # create Threshold_2 as the remaining thresholds
  dynamic_threshold = mock_data$P_Threshold
  
  # calculate differences in R2 between fixed and dynamic thresholds
  r2_differences = mock_data$R2 - max_r2_row$R2
  
  # build the data frame
  data.frame(
    Threshold_1 = fixed_threshold,
    Threshold_2 = dynamic_threshold,
    R2_Difference = r2_differences
  )
}


sensitivity_eur = sensitivity_analysis(mock_eur)
sensitivity_afr = sensitivity_analysis(mock_afr)
sensitivity_eas = sensitivity_analysis(mock_eas)

# statistical Tests
statistical_tests = function(y, genotypes) {
  # create a grouping factor based on the median of y (e.g., high vs low)
  group = factor(ifelse(y > median(y), "high", "low"))
  
  # Pearson correlation
  pearson_results = cor(y, genotypes[, 1], method = "pearson")
  
  # Spearman correlation
  spearman_results = cor(y, genotypes[, 1], method = "spearman")
  
  # Return results as a list (no trailing comma)
  list(
    Pearson = pearson_results,
    Spearman = spearman_results
  )
}



tests_eur = statistical_tests(phenotypes_eur$y, genotypes_eur)
tests_afr = statistical_tests(phenotypes_afr$y, genotypes_afr)
tests_eas = statistical_tests(phenotypes_eas$y, genotypes_eas)

# interpretation
interpret_results = function(tests, sensitivity) {
  cat("Statistical Test Results:\n")
  cat("Pearson Correlation:", tests$Pearson, "\n")
  cat("Spearman Correlation:", tests$Spearman, "\n")
  
  cat("\nSensitivity Analysis:\n")
  print(head(sensitivity))
}

interpret_results(tests_eur, sensitivity_eur)
interpret_results(tests_afr, sensitivity_afr)
interpret_results(tests_eas, sensitivity_eas)
