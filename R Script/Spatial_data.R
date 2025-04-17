# Spatial Data Analysis
# Used to check if the spatial distribution & dividing line makes sense
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(ggsignif)
library(ggpubr)
library(tidyverse)

spatial_data <- read.csv("Z:/Cembrowski Lab/MiaZ/BIOL449/Final_Ver/spatial_data.csv")
geom_data <- read.csv("Z:/Cembrowski Lab/MiaZ/BIOL449/Final_Ver/geometry_metadata.csv")

spatial_data <- spatial_data %>%
  mutate(side = case_when(
    side == "LEFT" ~ "Proximal",
    side == "RIGHT" ~ "Distal",
    TRUE ~ side  # other values remain unchanged
  ))

spatial_data <- spatial_data %>%
  filter(mouse_id != 1179205)

mouse_lookup <- data.frame(
  mouse_id = c("1142874", "1179201", "1179204", "1179202", 
               "1095199", "1095200", "1095202", "1095203"),
  mouse_label = c("01-6w-f", "04-10w-m", "05-10w-m", "06-6m-f",
                  "08-6m-m", "09-6m-m", "10-6m-m", "11-6m-m"),
  stringsAsFactors = FALSE
)

spatial_data$mouse_id <- as.character(spatial_data$mouse_id)

spatial_data <- spatial_data %>% 
  left_join(mouse_lookup, by = "mouse_id")

spatial_data <- spatial_data %>%
  arrange(mouse_id)

# 1. Simplified Cell Distribution Diagram in each SUB

# Merge using the unique identifier: mouse_id
data_merged <- merge(spatial_data, geom_data, by = "mouse_id")
cell_data_merged <- data_merged %>% filter(type == "cell")
ab_data_merged <- data_merged %>% filter(type == "plaque")

cell_side_colors <- c("Proximal" = "orange", "Distal" = "green4")
ab_side_colors <- c("Proximal" = "#FF82AB", "Distal" = "#5CACEE")

ggplot(cell_data_merged, aes(x = x_um, y = y_um, color = side)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = cell_side_colors) +
  labs(x = "X (um)", y = "Y (um)", title = "Spatial Distribution of Cells") +
  facet_wrap(~mouse_label, ncol = 4) +
  theme_minimal() +
  geom_abline(aes(slope = midline_slope, intercept = midline_intercept), 
              linetype = "solid", size = 1, color = "black") +
  geom_abline(aes(slope = divline_slope, intercept = divline_intercept), 
              linetype = "solid", size = 1, color = "red") +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 3)))+
  coord_fixed(
    ratio = 1, # ensure the scales of x and y are the same
    ylim = c(max(data_merged$y_um), 0) #(0, 0) in at top left corner
  )

ggplot(ab_data_merged, aes(x = x_um, y = y_um, color = side)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = ab_side_colors) +
  labs(x = "X (um)", y = "Y (um)", title = "Spatial Distribution of Plaques") +
  facet_wrap(~mouse_label, ncol = 4) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 3)))+
  geom_abline(aes(slope = midline_slope, intercept = midline_intercept), 
              linetype = "solid", size = 1, color = "black") +
  geom_abline(aes(slope = divline_slope, intercept = divline_intercept), 
              linetype = "solid", size = 1, color = "red") +
  theme(legend.position = "top") +
  coord_fixed(
    ratio = 1, # ensure the scales of x and y are the same
    ylim = c(max(data_merged$y_um), 0) #(0, 0) in at top left corner
  )
