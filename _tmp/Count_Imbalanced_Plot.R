library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)

# -----------------------------------------------------------------------------
# 1. RESHAPE AND PREPARE DATA
# -----------------------------------------------------------------------------
plot_df <- results_continuous %>%
  filter(!is.na(MIMOSA2_prob), !is.na(DiD_GLM_prob)) %>%
  pivot_longer(
    cols = c(MIMOSA2_prob, DiD_GLM_prob),
    names_to = "Method",
    values_to = "Probability"
  ) %>%
  mutate(
    Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD GLM"),
    Count_Level = factor(
      sub("_Count_.*", "", Cell_range),
      levels = c("1.00", ".90", ".75", ".50", ".25", ".10", ".050", ".025")
    ),
    Condition = ifelse(grepl("_Imb$", Cell_range), "Imbalanced", "Balanced")
  )

# -----------------------------------------------------------------------------
# 2. CALCULATE ROC COORDINATES PER GROUP
# -----------------------------------------------------------------------------
roc_coords <- plot_df %>%
  group_by(Count_Level, Condition, Method) %>%
  group_modify(~ {
    roc_obj <- pROC::roc(.x$Truth, .x$Probability, quiet = TRUE)
    data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities
    )
  }) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 3. PLOT 2x4 GRID OF ROC CURVES (2 ROWS, 4 COLUMNS)
# -----------------------------------------------------------------------------
ggplot(roc_coords, aes(x = FPR, y = TPR, color = Condition, linetype = Method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  geom_path(linewidth = 0.8, alpha = 0.85) +
  facet_wrap(~ Count_Level, nrow = 2, ncol = 4, labeller = label_both) +
  scale_color_manual(values = c("Balanced" = "#1F77B4", "Imbalanced" = "#D62728")) +
  scale_linetype_manual(values = c("MIMOSA2" = "solid", "DiD GLM" = "dashed")) +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Cell Balance",
    linetype = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )

#=================================================================================
library(dplyr)
library(ggplot2)
library(purrr)

# If not already installed:
# install.packages("pROC")
library(pROC)


# ==============================================================================
# 1. Prepare the data
# ==============================================================================

auc_data <- results_continuous %>%
  mutate(
    Count_Level = case_when(
      Cell_range == "1.00_Count_Ref"   ~ "1.00",
      grepl("^\\.90_Count", Cell_range)  ~ ".90",
      grepl("^\\.75_Count", Cell_range)  ~ ".75",
      grepl("^\\.50_Count", Cell_range)  ~ ".50",
      grepl("^\\.25_Count", Cell_range)  ~ ".25",
      grepl("^\\.10_Count", Cell_range)  ~ ".10",
      grepl("^\\.050_Count", Cell_range) ~ ".050",
      grepl("^\\.025_Count", Cell_range) ~ ".025",
      TRUE ~ NA_character_
    ),
    
    Balance = case_when(
      Cell_range == "1.00_Count_Ref" ~ "Balanced",
      grepl("_Bal$", Cell_range)     ~ "Balanced",
      grepl("_Imb$", Cell_range)     ~ "Imbalanced",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Count_Level), !is.na(Balance))


# ==============================================================================
# 2. Calculate AUC and analytic SE
#
#    We pool the observations within each Count_Level × Balance × Method
#    combination.
# ==============================================================================

calculate_auc <- function(dat, probability) {
  
  dat <- dat %>%
    filter(
      !is.na(.data[[probability]]),
      !is.na(Truth)
    )
  
  roc_obj <- roc(
    response  = dat$Truth,
    predictor = dat[[probability]],
    direction = "<",
    quiet     = TRUE
  )
  
  auc_value <- as.numeric(auc(roc_obj))
  
  # DeLong analytic variance
  auc_var <- as.numeric(var(roc_obj, method = "delong"))
  auc_se  <- sqrt(auc_var)
  
  data.frame(
    AUC = auc_value,
    AUC_SE = auc_se
  )
}


auc_results <- bind_rows(
  
  auc_data %>%
    group_by(Count_Level, Balance) %>%
    group_split() %>%
    map_dfr(function(dat) {
      
      out <- calculate_auc(dat, "MIMOSA2_prob")
      
      out %>%
        mutate(
          Count_Level = dat$Count_Level[1],
          Balance = dat$Balance[1],
          Method = "MIMOSA"
        )
    }),
  
  auc_data %>%
    group_by(Count_Level, Balance) %>%
    group_split() %>%
    map_dfr(function(dat) {
      
      out <- calculate_auc(dat, "DiD_GLM_prob")
      
      out %>%
        mutate(
          Count_Level = dat$Count_Level[1],
          Balance = dat$Balance[1],
          Method = "DiD"
        )
    })
)


# ==============================================================================
# 3. Get the 1.00 reference AUC for each method
# ==============================================================================

reference_auc <- auc_results %>%
  filter(Count_Level == "1.00") %>%
  select(
    Method,
    Ref_AUC = AUC,
    Ref_AUC_SE = AUC_SE
  )


# ==============================================================================
# 4. Calculate AUC odds ratio
#
#    OR = [AUC/(1-AUC)] / [Ref_AUC/(1-Ref_AUC)]
#
#    Delta-method SE on log(OR):
#
#    Var(logit(AUC)) ≈ Var(AUC) / [AUC^2 (1-AUC)^2]
#
#    Var(log OR) = Var(logit(AUC)) + Var(logit(Ref_AUC))
#
# ==============================================================================

or_results <- auc_results %>%
  left_join(reference_auc, by = "Method") %>%
  mutate(
    
    # Odds of AUC
    AUC_Odds = AUC / (1 - AUC),
    Ref_AUC_Odds = Ref_AUC / (1 - Ref_AUC),
    
    # Odds ratio
    AUC_OR = AUC_Odds / Ref_AUC_Odds,
    
    # SE of logit(AUC)
    SE_logit_AUC =
      AUC_SE / (AUC * (1 - AUC)),
    
    # SE of logit(reference AUC)
    SE_logit_Ref =
      Ref_AUC_SE / (Ref_AUC * (1 - Ref_AUC)),
    
    # SE of log(OR)
    SE_log_OR =
      sqrt(SE_logit_AUC^2 + SE_logit_Ref^2),
    
    # 95% CI on log scale, then exponentiate
    OR_Lower = exp(
      log(AUC_OR) - qnorm(0.975) * SE_log_OR
    ),
    
    OR_Upper = exp(
      log(AUC_OR) + qnorm(0.975) * SE_log_OR
    )
  )


# Remove the reference category from the plot
plot_data <- or_results %>%
  filter(Count_Level != "1.00") %>%
  mutate(
    Count_Level = factor(
      Count_Level,
      levels = c(
        ".90", ".75", ".50", ".25",
        ".10", ".050", ".025"
      )
    ),
    
    Balance = factor(
      Balance,
      levels = c("Balanced", "Imbalanced")
    ),
    
    Method = factor(
      Method,
      levels = c("MIMOSA", "DiD")
    )
  )


# ==============================================================================
# 5. Deterministic horizontal jitter
#
#    Four series are deliberately offset so they don't sit directly on top
#    of each other.
# ==============================================================================

plot_data <- plot_data %>%
  mutate(
    x_numeric = as.numeric(Count_Level),
    
    x_jitter = case_when(
      Method == "MIMOSA" & Balance == "Balanced"   ~ x_numeric - 0.12,
      Method == "MIMOSA" & Balance == "Imbalanced" ~ x_numeric - 0.04,
      Method == "DiD"    & Balance == "Balanced"   ~ x_numeric + 0.04,
      Method == "DiD"    & Balance == "Imbalanced" ~ x_numeric + 0.12
    )
  )


# ==============================================================================
# 6. Plot
# ==============================================================================

ggplot(
  plot_data,
  aes(
    x = x_jitter,
    y = AUC_OR,
    group = interaction(Method, Balance),
    colour = Method,
    linetype = Balance
  )
) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  geom_errorbar(
    aes(
      ymin = OR_Lower,
      ymax = OR_Upper
    ),
    width = 0.025,
    linewidth = 0.6
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = seq_along(
      c(".90", ".75", ".50", ".25",
        ".10", ".050", ".025")
    ),
    labels = c(
      ".90", ".75", ".50", ".25",
      ".10", ".050", ".025"
    )
  ) +
  
  scale_y_log10() +
  
  scale_colour_manual(
    values = c(
      "MIMOSA" = "#0072B2",
      "DiD"    = "#D55E00"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "Balanced"   = "solid",
      "Imbalanced" = "dashed"
    )
  ) +
  
  labs(
    title = "AUC Odds Ratios Relative to the 1.00 Count Reference",
    subtitle = "95% analytic CIs calculated using the delta method on the log(AUC odds ratio) scale",
    x = "Count level",
    y = "AUC Odds Ratio\n(relative to 1.00 reference)",
    colour = "Method",
    linetype = "Count balance"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )