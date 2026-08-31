library(ggplot2)
library(plotROC)
library(dplyr)
library(tidyr)

setwd("C:/Github/comp-2026-mimosa2")
load('_simulations/Count_Imbalance.Rdata')

# -----------------------------------------------------------------------------
# 1. RESHAPE AND PREPARE DATA
# -----------------------------------------------------------------------------
ROC_data_prepared <- results_continuous %>%
  filter(!is.na(MIMOSA2_prob), !is.na(DiD_GLM_prob)) %>%
  pivot_longer(
    cols = c(MIMOSA2_prob, DiD_GLM_prob),
    names_to = "Method",
    values_to = "Score"
  ) %>%
  mutate(
    Method = ifelse(Method == "MIMOSA2_prob", "MIMOSA2", "DiD Baseline"),
    Count_Level = factor(
      sub("_Count_.*", "", Cell_range),
      levels = c("1.00", ".90", ".75", ".50", ".25", ".10", ".050", ".025"),
      labels = c("100% target", "90% target", "75% target", "50% target", 
                 "25% target", "10% target", "5% target", "2.5% target")
    ),
    Condition = factor(
      ifelse(grepl("_Imb$", Cell_range), "Imbalanced", "Balanced"),
      levels = c("Balanced", "Imbalanced")
    )
  )

# -----------------------------------------------------------------------------
# 2. PLOT WITH GEOM_ROC
# -----------------------------------------------------------------------------
cell_imbalance_plot <- ggplot(
  data = ROC_data_prepared,
  mapping = aes(
    d = Truth,
    m = Score,
    color = Condition,
    linetype = Method,
    group = interaction(Method, Condition)
  )
) +
  geom_roc(n.cuts = 0, size = 1, linealpha = 0.9) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey10",
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ Count_Level,
    nrow = 2,
    ncol = 4
  ) +
  scale_color_manual(
    values = c(
      "Balanced" = "orchid",
      "Imbalanced" = "darkgreen"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "MIMOSA2" = "solid",
      "DiD Baseline" = "dotted"
    )
  ) +
  guides(
    linetype = guide_legend(
      order = 1,
      title = "Model",
      override.aes = list(linewidth = 1.2, size = 2)
    ),
    color = guide_legend(
      order = 2,
      title = "Parent cell type",
      override.aes = list(linewidth = 1.2, size = 1.5)
    )
  ) +
  labs(
    title = "MIMOSA2 performance under parent cell count imbalance.",
    subtitle = "ROC analysis.",
    x = "1-Specificity",
    y = "Sensitivity"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

print(cell_imbalance_plot)
ggsave('_fig/cell_imbalance_plot.pdf',plot=cell_imbalance_plot,width=8,height=6)
#=================================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(pROC)
library(ggplot2)

# ==============================================================================
# 1. PREPARE AUC DATA
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
# 2. CALCULATE AUC AND ANALYTIC SE
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
          Method = "MIMOSA2"
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
          Method = "DiD Baseline"
        )
    })
)


# ==============================================================================
# 3. GET REFERENCE AUC (1.00 COUNT)
# ==============================================================================
reference_auc <- auc_results %>%
  filter(Count_Level == "1.00") %>%
  select(
    Method,
    Ref_AUC = AUC,
    Ref_AUC_SE = AUC_SE
  )


# ==============================================================================
# 4. CALCULATE AUC ODDS RATIOS AND DELTA-METHOD CIs
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
    SE_logit_AUC = AUC_SE / (AUC * (1 - AUC)),
    
    # SE of logit(reference AUC)
    SE_logit_Ref = Ref_AUC_SE / (Ref_AUC * (1 - Ref_AUC)),
    
    # SE of log(OR)
    SE_log_OR = sqrt(SE_logit_AUC^2 + SE_logit_Ref^2),
    
    # 95% CI on log scale, then exponentiate
    OR_Lower = exp(log(AUC_OR) - qnorm(0.975) * SE_log_OR),
    OR_Upper = exp(log(AUC_OR) + qnorm(0.975) * SE_log_OR)
  )


# Filter out reference category and format factors
plot_data <- or_results %>%
  filter(Count_Level != "1.00") %>%
  mutate(
    Count_Level = factor(
      Count_Level,
      levels = c(".90", ".75", ".50", ".25", ".10", ".050", ".025")
    ),
    Balance = factor(
      Balance,
      levels = c("Balanced", "Imbalanced")
    ),
    Method = factor(
      Method,
      levels = c("MIMOSA2", "DiD Baseline")
    )
  )


# ==============================================================================
# 5. DETERMINISTIC HORIZONTAL JITTER
# ==============================================================================
plot_data <- plot_data %>%
  mutate(
    x_numeric = as.numeric(Count_Level),
    x_jitter = case_when(
      Method == "MIMOSA2"      & Balance == "Balanced"   ~ x_numeric - 0.12,
      Method == "MIMOSA2"      & Balance == "Imbalanced" ~ x_numeric - 0.04,
      Method == "DiD Baseline" & Balance == "Balanced"   ~ x_numeric + 0.04,
      Method == "DiD Baseline" & Balance == "Imbalanced" ~ x_numeric + 0.12
    )
  )


# ==============================================================================
# 6. PLOT (FORMATTED MATCHING ROC PLOTS)
# ==============================================================================
count_labels <- c("90% target", "75% target", "50% target", "25% target", 
                  "10% target", "5% target", "2.5% target")

auc_or_plot <- ggplot(
  plot_data,
  aes(
    x = x_jitter,
    y = AUC_OR,
    group = interaction(Method, Balance),
    colour = Balance,
    linetype = Method
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey10",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 1,
    alpha = 0.85
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = seq_along(count_labels),
    labels = count_labels
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      "Balanced"   = "orchid",
      "Imbalanced" = "darkgreen"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "MIMOSA2"      = "solid",
      "DiD Baseline" = "dashed"
    )
  ) +
  guides(
    linetype = guide_legend(
      order = 1,  # Swapped to 1 (Model comes first)
      title = "Model",
      override.aes = list(linewidth = 1.2, size = 2)
    ),
    colour = guide_legend(
      order = 2,  # Swapped to 2 (Parent cell balance comes second)
      title = "Parent cell balance",
      override.aes = list(linewidth = 1.2, size = 1.5)
    )
  ) +
  labs(
    title = "MIMOSA2 performance under parent cell count imbalance.",
    subtitle = "AUC odds ratios relative to 100% target count reference.",
    x = "Target downsampled parent cell count",
    y = "AUC odds ratio",
    colour = "Parent cell balance",
    linetype = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1, "cm")
  )

print(auc_or_plot)
# ggsave('_fig/auc_or_plot.pdf',plot=auc_or_plot,width=8,height=6)
