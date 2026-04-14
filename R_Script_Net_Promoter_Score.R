#####
## R code for "Coverage and Precision of Net Promoter Score Confidence Intervals Across Sampling Distributions"
## Philip Turk, Jordan Cinderich, Emma McNeill
## 2026/04/15
#####

## Load required libraries
library(tidyverse)
library(boot)

## Set seed for reproducibility
set.seed(2025)

## Define populations and true NPS values
populations <- list(
  E  = list(p = c(0.5, 0.0, 0.5), true_nps = 0),
  LS = list(p = c(1/12, 1/12, 5/6), true_nps = 0.75),
  T  = list(p = c(0.25, 0.5, 0.25), true_nps = 0),
  U  = list(p = c(1/3, 1/3, 1/3), true_nps = 0)
)

sample_sizes <- c(20, 50, 100)
n_sim <- 100 ## Set to 100,000 for paper
B <- 1000 ## Can tweak

## Adjusted Wald pseudo-count weights
aw_weights <- list(
  AW_E  = c(3/2, 0,   3/2),
  AW_LS = c(1/4, 1/4, 5/2),
  AW_T  = c(3/4, 3/2, 3/4),
  AW_U  = c(1,   1,   1)
)

## Unbiased variance estimator for NPS-hat
var_unbiased <- function(x, n) {
  p_hat <- x / n
  nps_hat <- p_hat[3] - p_hat[1]
  var_nps <- (p_hat[3]*(1 - p_hat[3]) + p_hat[1]*(1 - p_hat[1]) + 2*p_hat[3]*p_hat[1]) / (n - 1)
  return(setNames(c(nps_hat, sqrt(var_nps)), c("est", "se")))
}

## Bootstrap t confidence interval
bootstrap_t_ci <- function(x, n, B = 1000, alpha = 0.05) {
  ## Original estimate and standard error using Theorem 5
  original <- var_unbiased(x, n)
  theta_hat <- original["est"]
  se_hat <- original["se"]
  
  ## Reconstruct original sample data as a vector of class labels
  categories <- rep(c("D", "N", "P"), times = x)
  
  ## Initialize vector to store t* statistics
  t_star <- numeric(B)
  
  for (b in 1:B) {
    ## Nonparametric bootstrap: resample with replacement
    sample_b <- sample(categories, size = n, replace = TRUE)
    
    ## Tabulate counts in the same order (D, N, P)
    x_star <- table(factor(sample_b, levels = c("D", "N", "P")))
    x_star <- as.numeric(x_star)
    
    ## Estimate NPS* and SE* from bootstrap sample
    star <- var_unbiased(x_star, n)
    est_star <- star["est"]
    se_star  <- star["se"]
    
    ## Compute t* (only if SE* > 0 to avoid division by zero)
    t_star[b] <- if (se_star > 0) {
      (est_star - theta_hat) / se_star
    } else {
      NA_real_
    }
  }
  
  ## Compute bootstrap quantiles (ignoring any NAs)
  t_quantiles <- quantile(t_star, probs = c(1 - alpha/2, alpha/2), na.rm = TRUE)
  
  ## Compute final bootstrap t interval using SE from original sample
  ci_lower <- theta_hat - t_quantiles[1] * se_hat
  ci_upper <- theta_hat - t_quantiles[2] * se_hat
  
  return(setNames(c(theta_hat, ci_lower, ci_upper), c("est", "lower", "upper")))
}

## Adjusted Wald CI using pseudo-counts
adjusted_ci <- function(x, w, n) {
  x_adj <- x + w
  n_adj <- n + sum(w)
  p_hat <- x_adj / n_adj
  nps_hat <- p_hat[3] - p_hat[1]
  var_hat <- (p_hat[3]*(1 - p_hat[3]) + p_hat[1]*(1 - p_hat[1]) + 2*p_hat[3]*p_hat[1]) / (n_adj - 1)
  se_hat <- sqrt(var_hat)
  ci <- c(nps_hat - 1.96 * se_hat, nps_hat + 1.96 * se_hat)
  return(setNames(c(nps_hat, ci[1], ci[2]), c("est", "lower", "upper")))
}

## Run simulations
results <- list()
sampling_dists <- list()

for (pop_name in names(populations)) {
  for (n in sample_sizes) {
    pop <- populations[[pop_name]]
    for (i in 1:n_sim) {
      x <- rmultinom(1, n, pop$p)
      x <- as.vector(x)
      names(x) <- c("D", "N", "P")
      true_nps <- pop$true_nps
      
      ## Save NPS-hat for distribution plot
      p_hat <- x / n
      nps_hat <- p_hat[3] - p_hat[1]
      sampling_dists[[length(sampling_dists) + 1]] <- tibble(
        Population = pop_name, SampleSize = n, NPS_hat = nps_hat
      )
      
      ## Wald
      wald <- var_unbiased(x, n)
      ci_wald <- c(wald["est"] - 1.96 * wald["se"], wald["est"] + 1.96 * wald["se"])
      
      ## Bootstrap-t
      bt <- tryCatch(bootstrap_t_ci(x, n), error = function(e) c(est = NA, lower = NA, upper = NA))
      
      ## Adjusted Wald
      aw_results <- map_dfr(aw_weights, ~{
        ci <- adjusted_ci(x, ., n)
        tibble(lower = ci["lower"], upper = ci["upper"], est = ci["est"])
      }, .id = "Method")
      
      ## Combine CI results
      temp <- tibble(
        Population = pop_name,
        SampleSize = n,
        Method = c("Wald", "Bootstrap-t", aw_results$Method),
        Est = c(wald["est"], bt["est"], aw_results$est),
        Lower = c(ci_wald[1], bt["lower"], aw_results$lower),
        Upper = c(ci_wald[2], bt["upper"], aw_results$upper),
        Coverage = c(
          between(true_nps, ci_wald[1], ci_wald[2]),
          between(true_nps, bt["lower"], bt["upper"]),
          mapply(between, true_nps, aw_results$lower, aw_results$upper)
        ),
        Width = c(
          diff(ci_wald),
          bt["upper"] - bt["lower"],
          aw_results$upper - aw_results$lower
        )
      )
      
      results[[length(results) + 1]] <- temp
    }
  }
}

df_results <- bind_rows(results)
df_sampling <- bind_rows(sampling_dists)

## Sampling distribution faceted plot
## Reorder factors
df_sampling$Population <- factor(df_sampling$Population, levels = c("E", "LS", "T", "U"))
df_sampling$SampleSize <- factor(df_sampling$SampleSize, levels = c(20, 50, 100))

## Verify unbiasedness
# df_sampling %>% group_by(Population, SampleSize) %>% summarize(mean(NPS_hat))

## Reference truth lines
truth_df <- expand.grid(
  Population = levels(df_sampling$Population),
  SampleSize = levels(df_sampling$SampleSize)
) %>%
  mutate(TrueNPS = sapply(Population, function(p) populations[[p]]$true_nps))

bins_fd <- nclass.Sturges(df_sampling$NPS_hat)

## Final plot
ggplot(df_sampling, aes(x = NPS_hat)) +
  geom_histogram(aes(y = after_stat(density)),
    bins = 20, fill = "grey70", color = "white", alpha = 0.75) +
  geom_density(adjust = 2.5) +
  geom_vline(data = truth_df, aes(xintercept = TrueNPS), 
             color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_grid(rows = vars(SampleSize), cols = vars(Population)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(
    title = "Sampling Distributions of NPS-hat",
    subtitle = "Rows: Sample Size, Columns: Population",
    x = "NPS-hat",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

## Save as PNG
ggsave("Plots/plot01.png", width = 10, height = 8, units = "in", dpi = 600)

## Summarize CI performance
summary_df <- df_results %>%
  mutate(TrueNPS = sapply(Population, function(p) populations[[p]]$true_nps)) %>%
  group_by(Population, SampleSize, Method) %>%
  summarize(Coverage = mean(Coverage, na.rm = TRUE),
            Avg_CI_Width = mean(Width, na.rm = TRUE),
            Bias = mean(Est - TrueNPS, na.rm = TRUE),
            .groups = "drop")
# write.table(summary_df, file = "summary_df.csv", row.names = FALSE, quote = FALSE, sep = ",")

## Line/point plot of coverage with custom shapes
ggplot(summary_df, aes(x = factor(SampleSize), y = Coverage,
                       group = Method,
                       color = Method,
                       linetype = Method,
                       shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~ Population) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 1) +
  labs(title = "Coverage of 95% Confidence Intervals for NPS",
       x = "Sample Size", y = "Coverage") +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(3, "cm"),
    legend.key.height = unit(0.6, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size = 4))) +
  scale_shape_manual(values = c(
    "AW_E"       = 16,  # filled circle
    "AW_LS"      = 17,  # filled triangle
    "AW_T"       = 15,  # filled square
    "AW_U"       = 18,  # filled diamond (much more visible than cross)
    "Bootstrap-t"= 0,   # open square
    "Wald"       = 8    # star
  ))

## Save as PNG
ggsave("Plots/plot02_lines.png", width = 10, height = 8, units = "in", dpi = 600)

## Line/point plot of width with custom shapes
ggplot(summary_df, aes(x = factor(SampleSize), y = Avg_CI_Width,
                       group = Method,
                       color = Method,
                       linetype = Method,
                       shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~ Population) +
  labs(title = "Average Width of 95% Confidence Intervals for NPS",
       x = "Sample Size", y = "Average CI Width") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(3, "cm"),
    legend.key.height = unit(0.6, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size = 4))) +
  scale_shape_manual(values = c(
    "AW_E"       = 16,  # filled circle
    "AW_LS"      = 17,  # filled triangle
    "AW_T"       = 15,  # filled square
    "AW_U"       = 18,  # filled diamond (much more visible than cross)
    "Bootstrap-t"= 0,   # open square
    "Wald"       = 8    # star
  ))

## Save as PNG
ggsave("Plots/plot03_lines.png", width = 10, height = 8, units = "in", dpi = 600)
