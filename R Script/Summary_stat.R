library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(ggsignif)
library(ggpubr)
library(tidyverse)

summary_data <- read.csv("Z:/Cembrowski Lab/MiaZ/BIOL449/Final_Ver/summary_stats.csv")
summary_data <- summary_data %>%
  filter(mouse_id != 1179205)

mouse_lookup <- data.frame(
  mouse_id = c("1142874", "1179201", "1179204", "1179202", 
               "1095199", "1095200", "1095202", "1095203"),
  mouse_label = c("01-6w-f", "04-10w-m", "05-10w-m", "06-6m-f",
                  "08-6m-m", "09-6m-m", "10-6m-m", "11-6m-m"),
  stringsAsFactors = FALSE
)

summary_data$mouse_id <- as.character(summary_data$mouse_id)

summary_data <- summary_data %>% 
  left_join(mouse_lookup, by = "mouse_id")

summary_data <- summary_data %>%
  mutate(side = case_when(
    side == "LEFT" ~ "Proximal",
    side == "RIGHT" ~ "Distal",
    side == "TOTAL" ~ "Total",
    TRUE ~ side  # other values remain unchanged
  ))
# summary_data$mouse_label <- paste(summary_data$mouse_id, summary_data$sex, summary_data$age, sep = "\n")
summary_data <- summary_data %>%
  arrange(mouse_label)
summary_data <- summary_data %>% mutate(age = factor(age, levels = c("6 weeks", "10 weeks", "6 months")))

sides <- c("Proximal", "Distal", "Total")
age_colors <-  c("6 weeks" = "lightsalmon", "10 weeks" = "palegreen2", "6 months" = "skyblue1")
age_colors_points <- c("6 weeks" = "coral", "10 weeks" = "green4", "6 months" = "royalblue3")
sex_colors <- c("Male" = "black", "Female" = "red3")
side_patterns <- c("Proximal" = "stripe", "Distal" = "crosshatch", "Total" = "none")
age_comparisons <- list(c("10 weeks", "6 weeks"), c("6 months", "6 weeks"), c("6 months", "10 weeks"))


# Plot the cell count bar graph
# each mouse has 3 bars representing the cell count in Proximal and Distal sides and the sum.
create_plot_cell_count <- function(data) {
  ggplot(data) +
    geom_bar_pattern(
      aes(x = mouse_label, 
          y = (n_cells / sub_area), 
          fill = age,          # color represents age
          pattern = side,      # pattern represents side
          color = sex
      ),
      stat = "identity",
      position = position_dodge(width = 0.9),
      pattern_angle = 45,
      pattern_density = 0.2,
      pattern_spacing = 0.03,
      pattern_key_scale_factor = 0.8,
      width = 0.7,
      linewidth = 1 # make border (representing sex) thicker
    ) +
    
    geom_text(
      aes(x = mouse_label, 
          y = (n_cells / sub_area), 
          label = n_cells,
          group = side),
      position = position_dodge(width = 0.9),
      vjust = -0.5,
      size = 3.5
    ) +
    
    theme_minimal() +
    labs(x = "Mouse ID and Age",
         y = "Cell Count Normalized by Subiculum Area (μm²)",
         title = "Cell Count Normalized by Subiculum Area",
         fill = "Age",
         pattern = "Side",
         color = "Sex"
         ) +
    
    scale_fill_manual(values = age_colors) +
    scale_pattern_manual(values = side_patterns) +
    scale_color_manual(values = sex_colors) +
    
    guides(
      fill = guide_legend(override.aes = list(pattern = "none")),  
      color = guide_legend(override.aes = list(pattern = "none", fill = "grey80")),
      pattern = guide_legend(override.aes = list(fill = "grey80"))
    ) +
        # Legends
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.5, "cm"))
}

plot_cell_count <- create_plot_cell_count(summary_data)
print(plot_cell_count)


mouse_colors <- c(
  # 6 weeks（1）
  "01-6w-f" = "#ffe976",
  
  # 10 weeks（2）
  "04-10w-m" = "#87CEEB",
  "05-10w-m" = "#4682B4",
  
  # 6 months
  "06-6m-f" = "#6BAF8F",
  "08-6m-m"   = "#7EB09B",
  "09-6m-m"   = "#8FB28A",
  "10-6m-m"   = "#A0B577",
  "11-6m-m"   = "#B1B864"
)

library(effectsize)
library(broom)

# Compare the cell count in each side, grouped by Age
create_cell_count_side_plot <- function(data) {
  for (current_side in sides) {
    
    subset_data <- summary_data %>% filter(side == current_side)
    
    # 1. ANOVA
    aov_model <- aov(n_cells / sub_area ~ age, data = subset_data)
    eta_sq <- eta_squared(aov_model)
    anova_results <- summary(aov_model)[[1]]
    final_results <- data.frame(
      Effect = "Age",
      DF_num = anova_results$Df[1],
      DF_den = anova_results$Df[2],
      F_value = anova_results$`F value`[1],
      p_value = anova_results$`Pr(>F)`[1],
      eta_sq = eta_sq$Eta2
    )
    
    
    # 2. Tukey
    tukey_result <- TukeyHSD(aov_model)
    tukey_df <- as.data.frame(tukey_result$age)
    print(current_side)
    print(final_results)
    print(tukey_df)
    
    
    # 3. significance marks
    annotations <- sapply(age_comparisons, function(pair) {
      row_idx <- which(rownames(tukey_df) == paste(pair, collapse="-"))
      if (length(row_idx) > 0) {
        p_val <- tukey_df$`p adj`[row_idx]
        if (p_val < 0.001) {
          "***"
        } else if (p_val < 0.01) {
          "**"
        } else if (p_val < 0.05) {
          "*"
        } else {
          "NS"
        }
      } else {
        "N/A"
      }
    })
    # 4. plot
    ggp <- ggplot(subset_data, aes(x = age, y = n_cells / sub_area)) +
      
      labs(x = "Age Group",
           y = "Cell Count / Subiculum Area (μm²)",
           title = paste("Cell Count (", current_side, ") Normalized by Subiculum Area")) +
      
      scale_pattern_manual(values = side_patterns) +
      scale_color_manual(values = sex_colors) +
      scale_fill_manual(values = mouse_colors) +
      
      guides(
        fill = guide_legend(override.aes = list(pattern = "none")),
        color = guide_legend(override.aes = list(pattern = "none", fill = "grey80")),
        pattern = guide_legend(override.aes = list(fill = "grey80"))
      ) +
      
      geom_bar_pattern(
        stat = "identity",
        position = position_dodge(width = 0.8),
        aes(fill = mouse_label, pattern = side),
        width = 0.7,
        pattern_angle = 45,
        pattern_density = 0.2,
        pattern_spacing = 0.03,
        pattern_key_scale_factor = 0.8,
        linewidth = 1
      ) +
      
      geom_signif(comparisons = age_comparisons,
                  annotations = annotations,
                  map_signif_level = TRUE,
                  step_increase = 0.15,
                  textsize = 4,
                  tip_length = 0.01) +
      ylim(0, 0.00125) +
      theme_minimal() +
      theme(legend.position = "top")
    
    print(ggp)
  }
}

create_cell_count_side_plot(summary_data)



# Compare LEFT (Proximal) vs RIGHT (Distal) CELL COUNT in 6 months ONLY (paried)

proximal_data_cell_6mon <- summary_data %>% 
  filter(side == "Proximal", age == "6 months") %>% 
  select(mouse_label, age, n_cells, sub_area) %>% 
  mutate(proximal_cell_count = n_cells / sub_area)

distal_data_cell_6mon <- summary_data %>% 
  filter(side == "Distal", age == "6 months") %>% 
  select(mouse_label, age, n_cells, sub_area) %>% 
  mutate(distal_cell_count = n_cells / sub_area)

paired_data_cell_6mon <- proximal_data_cell_6mon %>% 
  inner_join(distal_data_cell_6mon, 
             by = c("mouse_label", "age"),
             suffix = c("_proximal", "_distal"))

# normality check
# shapiro_test_left_cell_6mon <- shapiro.test(paired_data_cell_6mon$proximal_cell_count)
# # normal distribution, so it will perform paried t test
# shapiro_test_right_cell_6mon <- shapiro.test(paired_data_cell_6mon$distal_cell_count)
# 
# if (shapiro_test_left_cell_6mon$p.value > 0.05 & shapiro_test_right_cell_6mon$p.value > 0.05) {
#   test_result <- t.test(paired_data_cell_6mon$proximal_cell_count, paired_data_cell_6mon$distal_cell_count, paired = TRUE)
# } else {
#   test_result <- wilcox.test(paired_data_cell_6mon$proximal_cell_count, paired_data_cell_6mon$distal_cell_count, paired = TRUE)
# }

test_result <- t.test(paired_data_cell_6mon$proximal_cell_count, paired_data_cell_6mon$distal_cell_count, paired = TRUE)

print(test_result)

long_data_cell_6mon <- summary_data %>% 
  filter(side %in% c("Proximal", "Distal"), age == "6 months") %>%  # 添加年龄过滤
  mutate(side = factor(side, levels = c("Proximal", "Distal")))

ggplot(long_data_cell_6mon, aes(x = side, y = n_cells / sub_area, fill = side)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(aes(group = mouse_label, color = factor(age)),  # 用年龄分颜色
             position = position_dodge(width = 0.8), 
             size = 3) +  
  geom_text(
    aes(
      label = mouse_label, 
      group = mouse_label
    ), 
    position = position_dodge(width = 0.8),
    vjust = -1,
    hjust = 0.5,
    size = 5
  )+
  geom_line(aes(group = mouse_label), 
            position = position_dodge(width = 0.8), 
            color = "gray50", alpha = 0.5) +  
  labs(
    x = "Side",
    y = "Cell Count / Subiculum Area (μm2)",
    title = "Mean Normalized Cell Count in Proximal and Distal Subiculum",
    color = "Age"
  ) +
  scale_fill_manual(values = c("Proximal" = "lightsalmon", "Distal" = "palegreen")) +
  scale_color_manual(values = age_colors_points)+
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  geom_signif(
    comparisons = list(c("Proximal", "Distal")),
    annotations = sapply(test_result$p.value, function(p) {
      if (is.na(p)) return("N/A")
      if (p < 0.001) "***"
      else if (p < 0.01) "**"
      else if (p < 0.05) "*"
      else "NS"
    }),
    y_position = max(long_data_cell_6mon$n_cells / long_data_cell_6mon$sub_area) * 1.1,
    tip_length = 0.01
  )



# For Amyloid Beta below:

# Plot the plaque area bar graph
# each mouse has 3 bars representing the cell count in left
create_plot_ab_area <- function(data) {
  ggplot(data) +
    geom_bar_pattern(
      aes(x = mouse_label, 
          y = plaque_area_pct, 
          fill = age,          # color represents age
          pattern = side,      # pattern represents side
          color = sex          # color represents sex
      ),
      stat = "identity",
      position = position_dodge(width = 0.9),  # 并排显示柱子
      pattern_angle = 45,
      pattern_density = 0.2,
      pattern_spacing = 0.03,
      pattern_key_scale_factor = 0.8,
      width = 0.7,
      linewidth = 1  # make border (representing sex) thicker
    ) +
    
    geom_text(
      aes(x = mouse_label, 
          y = plaque_area_pct, 
          label = sprintf("%.2f%%", plaque_area_pct * 100),
          group = side),
      position = position_dodge(width = 0.9),
      vjust = -0.5,
      size = 3.5
    ) +
    
    theme_minimal() +
    labs(x = "Mouse ID and Age",
         y = "Total Plaque Area Normalized by Subiculum Area",
         title = "Percentage of Total Plaque in Subiculum",
         fill = "Age",
         pattern = "Side",
         color = "Sex"
    ) +
    
    scale_fill_manual(values = age_colors) +
    scale_pattern_manual(values = side_patterns) +
    scale_color_manual(values = sex_colors) +
    
    guides(
      fill = guide_legend(override.aes = list(pattern = "none")),  
      color = guide_legend(override.aes = list(pattern = "none", fill = "grey80")),
      pattern = guide_legend(override.aes = list(fill = "grey80"))
    ) +
    
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.5, "cm")
    )
}

plot_ab_area <- create_plot_ab_area(summary_data)
print(plot_ab_area)

# Compare the cell count in each side, grouped by Age
create_plaque_area_side_plot <- function(data) {
  for (current_side in sides) {
    
    subset_data <- summary_data %>% filter(side == current_side)

    # 1. ANOVA
    aov_model <- aov(plaque_area_pct ~ age, data = subset_data)
    eta_sq <- eta_squared(aov_model)
    anova_results <- summary(aov_model)[[1]]
    final_results <- data.frame(
      Effect = "Age",
      DF_num = anova_results$Df[1],
      DF_den = anova_results$Df[2],
      F_value = anova_results$`F value`[1],
      p_value = anova_results$`Pr(>F)`[1],
      eta_sq = eta_sq$Eta2
    )
    
    
    # 2. Tukey
    tukey_result <- TukeyHSD(aov_model)
    tukey_df <- as.data.frame(tukey_result$age)
    print(current_side)
    print(final_results)
    print(tukey_df)
    
    
    # 3. significance marks
    annotations <- sapply(age_comparisons, function(pair) {
      row_idx <- which(rownames(tukey_df) == paste(pair, collapse="-"))
      if (length(row_idx) > 0) {
        p_val <- tukey_df$`p adj`[row_idx]
        if (p_val < 0.001) {
          "***"
        } else if (p_val < 0.01) {
          "**"
        } else if (p_val < 0.05) {
          "*"
        } else {
          "NS"
        }
      } else {
        "N/A"
      }
    })
    
    # 4. plot
    ggp <- ggplot(subset_data, aes(x = age, y = plaque_area_pct)) +
      
      labs(x = "Age Group",
           y = "Percentage of Plaque Area in Subiculum",
           title = paste("Plaque Area (", current_side, ") Normalized by Subiculum Area")) +
      
      scale_pattern_manual(values = side_patterns) +
      scale_color_manual(values = sex_colors) +
      scale_fill_manual(values = mouse_colors) +
      
      guides(
        fill = guide_legend(override.aes = list(pattern = "none")),
        color = guide_legend(override.aes = list(pattern = "none", fill = "grey80")),
        pattern = guide_legend(override.aes = list(fill = "grey80"))
      ) +
      
      geom_bar_pattern(
        stat = "identity",
        position = position_dodge(width = 0.8),
        aes(fill = mouse_label, pattern = side),
        width = 0.7,
        pattern_angle = 45,
        pattern_density = 0.2,
        pattern_spacing = 0.03,
        pattern_key_scale_factor = 0.8,
        linewidth = 1
      ) +
      
      geom_signif(comparisons = age_comparisons,
                  annotations = annotations,
                  map_signif_level = TRUE,
                  step_increase = 0.15,
                  textsize = 4,
                  tip_length = 0.01) +
      ylim(0, 0.3) +
      
      theme_minimal() +
      theme(legend.position = "top")
    
    print(ggp)
  }
}
create_plaque_area_side_plot(summary_data)


# Compare LEFT (Proximal) vs RIGHT (Distal) Plaque in 6 months ONLY (paried)

proximal_data_6mon <- summary_data %>% 
  filter(side == "Proximal", age == "6 months") %>%
  select(mouse_label, age, plaque_area_pct) %>% 
  rename(left_pct = plaque_area_pct)

distal_data_6mon <- summary_data %>% 
  filter(side == "Distal", age == "6 months") %>%
  select(mouse_label, age, plaque_area_pct) %>% 
  rename(right_pct = plaque_area_pct)

# Pair up
paired_data_6mon <- proximal_data_6mon %>% 
  inner_join(distal_data_6mon, by = c("mouse_label", "age"))

shapiro_test_left_6mon <- shapiro.test(paired_data_6mon$left_pct)
# data:  paired_data_6mon$left_pct
# W = 0.90417, p-value = 0.3992
# normal distribution, so it will perform paried t test
# shapiro_test_right_6mon <- shapiro.test(paired_data_6mon$right_pct)
# 
# if (shapiro_test_left_6mon$p.value > 0.05 & shapiro_test_right_6mon$p.value > 0.05) {
#   test_result <- t.test(paired_data_6mon$left_pct, paired_data_6mon$right_pct, paired = TRUE)
# } else {
#   test_result <- wilcox.test(paired_data_6mon$left_pct, paired_data_6mon$right_pct, paired = TRUE)
# }

test_result <- t.test(paired_data_6mon$left_pct, paired_data_6mon$right_pct, paired = TRUE)
print(test_result)

long_data_6mon <- summary_data %>% 
  filter(side %in% c("Proximal", "Distal"), age == "6 months") %>%
  mutate(side = factor(side, levels = c("Proximal", "Distal")))

ggplot(long_data_6mon, aes(x = side, y = plaque_area_pct, fill = side)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(aes(group = mouse_label, color = factor(age)),
             position = position_dodge(width = 0.8), 
             size = 3) +
  geom_text(
    aes(
      label = mouse_label, 
      group = mouse_label
    ), 
    position = position_dodge(width = 0.8),
    vjust = -1,
    hjust = 0.5,
    size = 5
  )+
  geom_line(aes(group = mouse_label), 
            position = position_dodge(width = 0.8), 
            color = "gray50", alpha = 0.5) +  
  labs(
    x = "Side",
    y = "Plaque Area Percentage",
    title = "Mean Plaque Area in Proximal and Distal Subiculum",
    color = "Age"
  ) +
  scale_fill_manual(values = c("Proximal" = "pink", "Distal" = "lightblue")) +
  scale_color_manual(values = age_colors_points)+
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  geom_signif(
    comparisons = list(c("Proximal", "Distal")),
    annotations = sapply(test_result$p.value, function(p) {
      if (is.na(p)) return("N/A")
      if (p < 0.001) "***"
      else if (p < 0.01) "**"
      else if (p < 0.05) "*"
      else "NS"
    }),
    y_position = max(long_data_6mon$plaque_area_pct) * 1.1,
    tip_length = 0.01
  )

# AB & Cell Correlation Analysis
library(ggrepel)

plot_normalized_regression <- function(side_choice) {
  
  temp_data <- summary_data %>% 
    filter(side == side_choice) %>% 
    mutate(age = factor(age, levels = c("6 weeks", "10 weeks", "6 months")))

  temp_data$normalized_n_cells_sub_area <- scale(temp_data$n_cells / temp_data$sub_area)
  temp_data$normalized_plaque_area_pct <- scale(temp_data$plaque_area_pct)

  model <- lm(normalized_n_cells_sub_area ~ normalized_plaque_area_pct, data = temp_data)
  print(summary(model))
  sd_plaque <- sd(temp_data$plaque_area_pct)
  sd_neuronal <- sd(temp_data$n_cells / temp_data$sub_area)
  print(paste("SD for AB:", sd_plaque))
  print(paste("SD for Neuron:", sd_neuronal))
  
  r_squared <- summary(model)$r.squared
  p_value <- summary(model)$coefficients[2, 4]
  y_max <- max(temp_data$normalized_n_cells_sub_area)
  
  ggplot(temp_data, aes(x = normalized_plaque_area_pct, y = normalized_n_cells_sub_area)) +
    geom_point(aes(color = age, shape = sex), size = 4, alpha = 0.7) +  # Color by Age and shape by Sex
    geom_text_repel(
      aes(label = mouse_label),
      size = 5,
      box.padding = 0.5,  # Space between label and point
      direction = "both"
    )+
    geom_smooth(method = "lm", se = TRUE,
                linetype = "solid", size = 1, 
                color = "black",
                fill = "grey60",
                alpha = 0.2) +
    labs(x = "Normalized Plaque Area in Subiculum",
         y = "Normalized Cell Count / Subiculum Area",
         title = paste("Normalized Neuron Count vs. Normalized Plaque Area (", side_choice, ")"),
         color = "Age",
         shape = "Sex") +
    theme_minimal() +
    theme(
      # Customize legend text and titles:
      legend.text = element_text(size = 14),    # Legend item labels (e.g., "6 weeks", "Female")
      legend.title = element_text(size = 16, face = "bold")  # Legend titles ("Age", "Sex")
    )+
    expand_limits(y = y_max * 1.1) +
    scale_color_manual(values = age_colors_points,
                       breaks = c("6 weeks", "10 weeks", "6 months")) +  # Explicit order of age
    scale_shape_manual(values = c("Female" = 19, "Male" = 15)) +
    annotate("text", x = 0.5, y = y_max * 0.8, 
             label = paste("R² =", round(r_squared, 3)), size = 7, color = "black") +
    annotate("text", x = 0.5, y = y_max * 0.95, 
             label = paste("p =", format(p_value, digits = 2)), size = 7, color = "black")
}

plot_normalized_regression("Total")
plot_normalized_regression("Proximal")
plot_normalized_regression("Distal")