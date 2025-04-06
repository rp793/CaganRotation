# Code to create all three CSV files from the raw data

# 1. Create parental phased data CSV
write_parental_csv <- function() {
  data <- data.frame(
    Sample = c("GS-1_ERR3284981", "GS-1_ERR3284982", "MS-1_ERR11203029", "MS-2_ERR11203031", "MS-2_ERR11203032"),
    Total_DNMs = c(658, 640, 832, 938, 818),
    Father_Phased = c(467, 447, 590, 613, 520),
    Mother_Phased = c(191, 193, 242, 325, 298)
  )
  
  write.csv(data, "parental-phased-data.csv", row.names = FALSE)
  cat("Created parental-phased-data.csv\n")
}
# 2. Create mutation types data CSV
write_mutation_csv <- function() {
  data <- data.frame(
    Sample = c("GS-1_ERR3284981", "GS-1_ERR3284982", "MS-1_ERR11203029", "MS-2_ERR11203031", "MS-2_ERR11203032"),
    A_C = c(44, 41, 51, 40, 38),
    A_G = c(95, 69, 92, 113, 95),
    A_T = c(27, 31, 57, 59, 55),
    C_A = c(29, 42, 48, 67, 46),
    C_G = c(20, 38, 40, 46, 47),
    C_T = c(94, 103, 142, 148, 115),
    G_A = c(120, 124, 139, 137, 127),
    G_C = c(31, 21, 37, 51, 49),
    G_T = c(34, 37, 39, 54, 43),
    T_A = c(41, 33, 46, 45, 61),
    T_C = c(79, 79, 104, 120, 103),
    T_G = c(44, 22, 37, 58, 39)
  )
  
  write.csv(data, "mutation-types-data.csv", row.names = FALSE)
  cat("Created mutation-types-data.csv\n")
}

# 3. Create chromosome distribution data CSV
write_chromosome_csv <- function() {
  # Create data frame with chromosome numbers and mutation counts
  chrom_data <- data.frame(
    Chromosome = c(1:38, "X"),
    GS_1_ERR3284981 = c(31, 18, 17, 20, 10, 23, 13, 23, 12, 6, 12, 9, 8, 60, 20, 12, 
                        8, 11, 93, 4, 9, 3, 12, 11, 13, 52, 9, 3, 18, 4, 11, 21, 3, 
                        7, 7, 4, 3, 5, 53),
    GS_1_ERR3284982 = c(15, 10, 10, 8, 16, 28, 14, 18, 2, 10, 6, 14, 10, 44, 41, 34, 
                        19, 2, 64, 14, 24, 7, 16, 6, 0, 43, 5, 6, 11, 6, 1, 22, 22, 
                        7, 5, 3, 9, 8, 60),
    MS_1_ERR11203029 = c(11, 14, 15, 11, 5, 14, 4, 43, 2, 10, 12, 13, 6, 23, 39, 19, 
                         23, 11, 125, 15, 19, 11, 11, 20, 11, 29, 39, 14, 19, 9, 13, 68, 
                         4, 22, 3, 16, 1, 3, 105),
    MS_2_ERR11203031 = c(64, 8, 24, 11, 8, 33, 7, 56, 13, 16, 16, 11, 13, 61, 75, 12, 
                         12, 8, 109, 7, 19, 16, 27, 12, 13, 29, 10, 19, 4, 8, 16, 43, 
                         12, 18, 15, 6, 2, 6, 99),
    MS_2_ERR11203032 = c(14, 15, 22, 10, 12, 29, 13, 43, 16, 18, 21, 4, 17, 61, 23, 16, 
                         22, 10, 85, 7, 10, 21, 10, 10, 6, 45, 22, 5, 8, 11, 7, 30, 13, 
                         13, 15, 13, 5, 1, 115)
  )
  
  # Rename columns to match the format in your original data
  names(chrom_data) <- c("Chromosome", "GS-1_ERR3284981", "GS-1_ERR3284982", "MS-1_ERR11203029", 
                         "MS-2_ERR11203031", "MS-2_ERR11203032")
  
  write.csv(chrom_data, "chromosome-data.csv", row.names = FALSE)
  cat("Created chromosome-data.csv\n")
}

# Run all the functions to create the CSV files
write_parental_csv()
write_mutation_csv()
write_chromosome_csv()

# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

# Read the data files
parental_data <- read.csv("parental-phased-data.csv")
mutation_data <- read.csv("mutation-types-data.csv")
chrom_data <- read.csv("chromosome-data.csv")

# 1. Parental Phased Variants Visualization
parental_data_long <- parental_data %>%
  pivot_longer(cols = c(Father_Phased, Mother_Phased),
               names_to = "Parent",
               values_to = "Count") %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / Total_DNMs)

# Create the plot for parental phased variants
ggplot(parental_data_long, aes(x = Sample, y = Proportion, fill = Parent)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set1", 
                    labels = c("Father_Phased" = "Paternal", "Mother_Phased" = "Maternal")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of Paternal vs Maternal Phased Variants",
       x = "Sample",
       y = "Proportion of Total DNMs",
       fill = "Parent")
ggsave("parental-phased-variants.png", width = 10, height = 6)

# 2. Mutation Types Visualization
mutation_data_long <- mutation_data %>%
  pivot_longer(cols = -Sample, 
               names_to = "Mutation_Type", 
               values_to = "Count") %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count),
         Proportion = Count / Total,
         Mutation_Type = gsub("_", ">", Mutation_Type))

# Create the plot for mutation types
ggplot(mutation_data_long, aes(x = Sample, y = Proportion, fill = Mutation_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Proportion of Mutation Types",
       x = "Sample",
       y = "Proportion",
       fill = "Mutation Type")
ggsave("mutation-types.png", width = 12, height = 6)

# 3. Chromosome Distribution Visualization
chrom_data_long <- chrom_data %>%
  pivot_longer(cols = -Chromosome, 
               names_to = "Sample", 
               values_to = "Count")

# Calculate proportions by sample
chrom_proportions <- chrom_data_long %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count))

# Create a heatmap for chromosome distribution
ggplot(chrom_proportions, aes(x = Sample, y = factor(Chromosome), fill = Proportion)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Mutations Across Chromosomes",
       x = "Sample",
       y = "Chromosome",
       fill = "Proportion")
ggsave("chromosome-heatmap.png", width = 10, height = 12)

# Create a chromosome count plot for the top 10 chromosomes
top_chroms <- chrom_data_long %>%
  group_by(Chromosome) %>%
  summarize(Total = sum(Count)) %>%
  arrange(desc(Total)) %>%
  head(10) %>%
  pull(Chromosome)

ggplot(chrom_data_long %>% filter(Chromosome %in% top_chroms), 
       aes(x = factor(Chromosome), y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Top 10 Chromosomes by Mutation Count",
       x = "Chromosome",
       y = "Count",
       fill = "Sample")
ggsave("top-chromosomes.png", width = 12, height = 6)

print("All visualizations created successfully!")

#adding barchart:
# 4. Mutation Types Bar Chart - absolute counts
# Reuse the mutation_data_long dataframe we created earlier

# Create the bar chart for mutation types (absolute counts)
ggplot(mutation_data_long, aes(x = Mutation_Type, y = Count, fill = Mutation_Type)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Sample, scales = "free_y") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none",
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Mutation Type Counts by Sample", 
    x = "Mutation Type", 
    y = "Count"
  )
ggsave("mutation-types-counts-by-sample.png", width = 12, height = 8)

# Create an aggregate bar chart across all samples
mutation_summary <- mutation_data_long %>%
  group_by(Mutation_Type) %>%
  summarize(Total_Count = sum(Count))

ggplot(mutation_summary, aes(x = reorder(Mutation_Type, -Total_Count), y = Total_Count, fill = Mutation_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Total Counts of Mutation Types Across All Samples",
    x = "Mutation Type",
    y = "Total Count"
  )
ggsave("mutation-types-total-counts.png", width = 10, height = 6)
