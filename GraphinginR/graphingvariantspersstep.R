# Load required libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Import the data
data_file <- "D:/PhD/rp793_Cagan/rotationresults/variantsateachstep.txt"
variant_data <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)

# Remove the Filtering_Known_Variants column as requested
variant_data <- variant_data %>% select(-Filtering_Known_Variants)

# Convert data to long format for easier plotting
variant_data_long <- variant_data %>%
  pivot_longer(
    cols = c(Called_SNPs_After_Hard_Filtering, De_Novo_Annotation, 
             Post_Filtering_Steps, DP20_VAF04.06_Filtering),
    names_to = "Filtering_Stage",
    values_to = "Variant_Count"
  )

# Factor the filtering stages for proper ordering in the plot
variant_data_long$Filtering_Stage <- factor(
  variant_data_long$Filtering_Stage,
  levels = c("Called_SNPs_After_Hard_Filtering", "De_Novo_Annotation", 
             "Post_Filtering_Steps", "DP20_VAF04.06_Filtering")
)

# Create more readable labels for the filtering stages
stage_labels <- c(
  "Called_SNPs_After_Hard_Filtering" = "Hard Filtered SNPs",
  "De_Novo_Annotation" = "De Novo Annotation",
  "Post_Filtering_Steps" = "Post Filtering",
  "DP20_VAF04.06_Filtering" = "DP20 VAF0.4-0.6"
)

# Create a grouped bar plot
ggplot(variant_data_long, aes(x = Family, y = Variant_Count, fill = Filtering_Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~Sample_ID, scales = "free_y") +
  scale_fill_brewer(palette = "Set1", labels = function(x) stage_labels[x]) +
  theme_minimal() +
  labs(
    title = "Variant Counts at Different Filtering Stages",
    x = "Family",
    y = "Number of Variants",
    fill = "Filtering Stage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Save the plot
ggsave("variant_filtering_stages.png", width = 12, height = 8, dpi = 300)

# Alternative visualization: Line plot showing the progression of filtering
ggplot(variant_data_long, aes(x = Filtering_Stage, y = Variant_Count, 
                              group = Sample_ID, color = Family)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_discrete(labels = function(x) stage_labels[x]) +
  theme_light() +
  labs(
    title = "Variant Reduction through Filtering Pipeline",
    x = "Filtering Stage",
    y = "Number of Variants",
    color = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_brewer(palette = "Dark2")

# Save the line plot
ggsave("variant_filtering_progression.png", width = 10, height = 6, dpi = 300)

# Create a summary table with mean variants per family at each stage
summary_by_family <- variant_data_long %>%
  group_by(Family, Filtering_Stage) %>%
  summarize(
    Mean_Variants = mean(Variant_Count, na.rm = TRUE),
    SD_Variants = sd(Variant_Count, na.rm = TRUE),
    .groups = "drop"
  )

# Print the summary table
print(summary_by_family)

# Create a plot for the summary by family
ggplot(summary_by_family, aes(x = Filtering_Stage, y = Mean_Variants, 
                              group = Family, color = Family)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_Variants - SD_Variants, 
                    ymax = Mean_Variants + SD_Variants), 
                width = 0.2) +
  scale_x_discrete(labels = function(x) stage_labels[x]) +
  theme_minimal() +
  labs(
    title = "Mean Variant Counts by Family through Filtering Stages",
    x = "Filtering Stage",
    y = "Mean Number of Variants",
    color = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the summary plot
ggsave("variant_family_summary.png", width = 10, height = 6, dpi = 300)

---
#reattempting:
# Load required libraries
  library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Import the data
data_file <- "D:/PhD/rp793_Cagan/rotationresults/variantsateachstep.txt"
variant_data <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)

# Remove the Filtering_Known_Variants column as requested
variant_data <- variant_data %>% select(-Filtering_Known_Variants)

# Remove any rows that might contain "variant_counts.csv (END)" in the Family column
variant_data <- variant_data %>% 
  filter(Family != "variant_counts.csv (END)")

# Convert data to long format for easier plotting
variant_data_long <- variant_data %>%
  pivot_longer(
    cols = c(Called_SNPs_After_Hard_Filtering, De_Novo_Annotation, 
             Post_Filtering_Steps, DP20_VAF04.06_Filtering),
    names_to = "Filtering_Stage",
    values_to = "Variant_Count"
  )

# Factor the filtering stages for proper ordering in the plot
variant_data_long$Filtering_Stage <- factor(
  variant_data_long$Filtering_Stage,
  levels = c("Called_SNPs_After_Hard_Filtering", "De_Novo_Annotation", 
             "Post_Filtering_Steps", "DP20_VAF04.06_Filtering")
)

# Create more readable labels for the filtering stages
stage_labels <- c(
  "Called_SNPs_After_Hard_Filtering" = "Hard Filtered SNPs",
  "De_Novo_Annotation" = "De Novo Annotation",
  "Post_Filtering_Steps" = "Post Filtering",
  "DP20_VAF04.06_Filtering" = "DP20 VAF0.4-0.6"
)

# Line plot showing the progression of filtering with adjusted y-axis
ggplot(variant_data_long, aes(x = Filtering_Stage, y = Variant_Count, 
                              group = Sample_ID, color = Family)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_discrete(labels = function(x) stage_labels[x]) +
  scale_y_continuous(
    trans = "log10",  # Use log scale to better see small values
    breaks = c(100, 200, 500, 1000, 10000, 100000, 1000000, 4000000),
    labels = scales::comma
  ) +
  annotation_logticks(sides = "l") +  # Add log ticks to the left side
  theme_light() +
  labs(
    title = "Variant Reduction through Filtering Pipeline",
    x = "Filtering Stage",
    y = "Number of Variants (log scale)",
    color = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2")

# Save the updated line plot
ggsave("variant_filtering_progression_updated.png", width = 10, height = 6, dpi = 300)

# Alternative visualization with linear scale but focus on lower range
ggplot(variant_data_long, aes(x = Filtering_Stage, y = Variant_Count, 
                              group = Sample_ID, color = Family)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_discrete(labels = function(x) stage_labels[x]) +
  scale_y_continuous(
    limits = c(0, 4500000),  # Set upper limit
    breaks = c(0, 200, 500, 1000, 10000, 100000, 1000000, 2000000, 3000000, 4000000),
    labels = scales::comma
  ) +
  theme_light() +
  labs(
    title = "Variant Reduction through Filtering Pipeline",
    x = "Filtering Stage",
    y = "Number of Variants",
    color = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_brewer(palette = "Dark2")

# Save the alternative visualization
ggsave("variant_filtering_progression_linear.png", width = 10, height = 6, dpi = 300)

# Create closeup view of just the final two filtering stages
final_stages_data <- variant_data_long %>%
  filter(Filtering_Stage %in% c("Post_Filtering_Steps", "DP20_VAF04.06_Filtering"))

ggplot(final_stages_data, aes(x = Filtering_Stage, y = Variant_Count, 
                              group = Sample_ID, color = Family)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = Variant_Count), vjust = -0.8, size = 3) +
  scale_x_discrete(labels = function(x) stage_labels[x]) +
  theme_light() +
  labs(
    title = "Final Filtering Stages - Close-up View",
    x = "Filtering Stage",
    y = "Number of Variants",
    color = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_brewer(palette = "Dark2")

# Save the closeup view
ggsave("variant_filtering_final_stages.png", width = 10, height = 6, dpi = 300)
