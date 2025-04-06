# Step 1: Create a CSV file from the hiConfDeNovo data

# Create data frame with the values from your file
hiconf_data <- data.frame(
  Family = c("GS-1", "GS-1", "MS-2", "MS-2", "MS-1"),
  Sample_ID = c("ERR3284981", "ERR3284982", "ERR11203031", "ERR11203032", "ERR11203029"),
  HiConfDeNovo_Count = c(23083, 22161, 22715, 20901, 19647)
)

# Save as CSV file
write.csv(hiconf_data, "hiconf_denovo_data.csv", row.names = FALSE)
cat("CSV file created: hiconf_denovo_data.csv\n")

# Step 2: Visualize the data

# Load libraries
library(ggplot2)

# Create a bar plot showing variant counts by sample
p1 <- ggplot(hiconf_data, aes(x = Sample_ID, y = HiConfDeNovo_Count, fill = Family)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = HiConfDeNovo_Count), vjust = -0.5, size = 3.5) +
  scale_fill_brewer(palette = "Set1") +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "High Confidence De Novo Variants by Sample",
    x = "Sample ID",
    y = "Number of Variants"
  )

# Save the plot
ggsave("hiconf_by_sample.png", p1, width = 9, height = 6)
cat("Plot saved: hiconf_by_sample.png\n")

# Create a bar plot grouped by family
p2 <- ggplot(hiconf_data, aes(x = Family, y = HiConfDeNovo_Count, fill = Sample_ID)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(
    aes(label = HiConfDeNovo_Count),
    position = position_dodge(0.8),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_light() +
  labs(
    title = "High Confidence De Novo Variants by Family",
    x = "Family",
    y = "Number of Variants",
    fill = "Sample ID"
  )

# Save the family plot
ggsave("hiconf_by_family.png", p2, width = 9, height = 6)
cat("Plot saved: hiconf_by_family.png\n")
