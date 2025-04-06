#analyse_denovo_updated
#for all samples no extra filtering:
library(ggplot2)
library(RColorBrewer)

# Load original data
all_Samples_snv_dnms <- read.table("D:/PhD/rp793_Cagan/rotationresults/all_denovo_snvs.txt", 
                                   header=FALSE, sep="\t", 
                                   col.names=c("sample", "chrom", "pos", "ref", "alt", "t_alt_count", "t_depth"))

# Load metadata
metadata <- read.delim("D:/PhD/rp793_Cagan/rotationresults/metadata.txt", 
                       header=TRUE, sep="\t")

# Merge with metadata, keeping only necessary columns
all_Samples_snv_dnms <- merge.data.frame(all_Samples_snv_dnms, 
                                         metadata[, c("Sample", "FamilyID", "FatherID", "MotherID", "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE")], 
                                         by.x="sample", by.y="Sample")

# Calculate DNM count per sample
for (i in unique(all_Samples_snv_dnms$sample)) {
  all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == i), "SBS_DNMs_Count"] <- 
    nrow(all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == i), ])
}

# Calculate VAF and coverage metrics
for (sample in unique(all_Samples_snv_dnms$sample)) {
  all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "sample.median.vaf"] <- 
    median(all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "t_alt_count"] /
             all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "t_depth"])
  all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "sample.mean.vaf"] <- 
    mean(all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "t_alt_count"] /
           all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "t_depth"])
  all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "Mean_Cov_variants"] <- 
    median(all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == sample), "t_depth"])
}

# Subset per-sample summary
all_Samples_snv_dnms_per_sample <- unique(all_Samples_snv_dnms[, c("sample", "FamilyID", "FatherID", "MotherID", 
                                                                   "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE", 
                                                                   "SBS_DNMs_Count", "sample.median.vaf", 
                                                                   "sample.mean.vaf", "Mean_Cov_variants")])

# Write table
write.table(all_Samples_snv_dnms_per_sample, 
            file="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_per_offspring.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

# Define colors based on FamilyID
nb.cols <- length(unique(all_Samples_snv_dnms_per_sample$FamilyID))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
names(mycolors) <- unique(all_Samples_snv_dnms_per_sample$FamilyID)

# Define colors for samples (for VAF plot)
nb.cols_samples <- length(unique(all_Samples_snv_dnms$sample))
sample_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols_samples)
names(sample_colors) <- unique(all_Samples_snv_dnms$sample)

# Plot 1: VAF Distribution (as histogram, color by sample, facet by FamilyID)
vaf_plot <- ggplot(all_Samples_snv_dnms, aes(x=t_alt_count/t_depth)) + 
  geom_histogram(aes(fill=sample), alpha=0.4, bins=30, position="identity") + 
  facet_wrap(vars(FamilyID), ncol=2, scales = "free_x") +  # Force each facet to have its own axis
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=14, margin=margin(b=10, t=10)),
        axis.text.y=element_text(size=11), 
        axis.text.x=element_text(angle=0, size=11, vjust=0.5),  # Simplified axis text positioning
        axis.title.x=element_text(size=12, margin=margin(t=30)),  # Increased top margin
        axis.title.y=element_text(size=14),
        plot.margin=unit(c(1, 1, 5, 1), "lines"),  # Significantly increased bottom margin
        plot.subtitle=element_text(hjust=0.5, size=14, margin=margin(t=10, b=10)),
        panel.spacing=unit(3, "lines"),  # Increased panel spacing
        # Remove custom axis line settings temporarily
        # axis.line.x = element_line(color="black", linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray95"),
        plot.title=element_text(hjust=0.5),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) + 
  labs(x="Variant Allele Frequency (VAF)", y="Variant Count", subtitle="Variant Allele Frequency (VAF) Distribution at De Novo Sites") +
  scale_fill_manual(name="Sample", values=sample_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), position = "bottom") +  # Force x-axis position
  coord_cartesian(clip = "off", expand = TRUE)  # Ensure no clipping

# Try rendering the plot with a specific aspect ratio
ggsave("vaf_plot.png", vaf_plot, width=12, height=10, dpi=300)

# Plot 2: DNM Count by Sample
dnm_by_sample <- ggplot(all_Samples_snv_dnms_per_sample, 
                        aes(x=sample, y=SBS_DNMs_Count, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  ylim(0, max(all_Samples_snv_dnms_per_sample$SBS_DNMs_Count) + 50) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Sample") + 
  ylab("De novo variants") + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text.x=element_text(angle=45, hjust=1, size=11),
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Plot 3: DNM Count vs. Mean Coverage
dnm_vs_coverage <- ggplot(all_Samples_snv_dnms_per_sample, 
                          aes(x=MEAN_COVERAGE, y=SBS_DNMs_Count, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, max(all_Samples_snv_dnms_per_sample$SBS_DNMs_Count) + 50) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Mean Coverage") + 
  ylab("De novo variants") + 
  geom_smooth(data=subset(all_Samples_snv_dnms_per_sample, 
                          FamilyID %in% names(which(table(all_Samples_snv_dnms_per_sample$FamilyID) >= 2))),
              aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.line.x=element_line(), 
        axis.line.y=element_line(),
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Plot 4: Distribution of Coverage for De Novo Variants (as histogram)
coverage_dist <- ggplot(all_Samples_snv_dnms, aes(x=t_depth)) + 
  geom_histogram(aes(fill=FamilyID), alpha=0.4, bins=30, position="identity") + 
  facet_wrap(vars(sample), ncol=2) + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=12, margin=margin(b=10, t=10)),
        axis.text=element_text(size=10), 
        axis.text.x=element_text(margin=margin(t=10)),
        axis.title.x=element_text(size=12, margin=margin(t=20)), 
        axis.title.y=element_text(size=12),
        plot.margin=unit(c(1, 1, 4, 1), "lines"),
        plot.title=element_text(hjust=0.5),
        panel.spacing=unit(2, "lines"),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray95")) +  
  xlab("Coverage (Depth) at De Novo Sites") + 
  ylab("Variant Count") + 
  scale_fill_manual(name="FamilyID", values=mycolors) + 
  xlim(0, 100) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  facet_wrap(vars(sample), ncol=2, scales = "free_x")  # Add free x scales

# Plot 5: Average Coverage vs. Average VAF
cov_vaf <- ggplot(all_Samples_snv_dnms_per_sample, 
                  aes(x=Mean_Cov_variants, y=sample.mean.vaf, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, 1) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Average Coverage at De Novo Sites") + 
  ylab("Average VAF") + 
  geom_smooth(data=subset(all_Samples_snv_dnms_per_sample, 
                          FamilyID %in% names(which(table(all_Samples_snv_dnms_per_sample$FamilyID) >= 2))),
              aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text=element_text(size=10), 
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Save all plots (increased width for VAF plot to accommodate legend)
ggsave(vaf_plot, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_VAF_distribution.png", 
       dpi=300, units="in", width=12, height=10)  # Increased width
ggsave(dnm_by_sample, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_by_Sample.png", 
       dpi=300, units="in", width=7, height=5)
ggsave(dnm_vs_coverage, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_vs_Coverage.png", 
       dpi=300, units="in", width=7, height=5)
ggsave(coverage_dist, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Coverage_Distribution_DeNovo.png", 
       dpi=300, units="in", width=10, height=10)
ggsave(cov_vaf, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Avg_Coverage_vs_Avg_VAF.png", 
       dpi=300, units="in", width=7, height=5)


#filtered file:
library(ggplot2)
library(RColorBrewer)

# Load filtered data
filtered_df_depth2 <- read.table("D:/PhD/rp793_Cagan/rotationresults/filtered_variants.txt", 
                                 header=FALSE, sep="\t", 
                                 col.names=c("sample", "chrom", "pos", "ref", "alt", "t_alt_count", "t_depth"))

# Load metadata
metadata <- read.delim("D:/PhD/rp793_Cagan/rotationresults/metadata.txt", 
                       header=TRUE, sep="\t")

# Merge with metadata, keeping only necessary columns
filtered_df_depth2 <- merge.data.frame(filtered_df_depth2, 
                                       metadata[, c("Sample", "FamilyID", "FatherID", "MotherID", "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE")], 
                                       by.x="sample", by.y="Sample")

# Calculate DNM count per sample
for (i in unique(filtered_df_depth2$sample)) {
  filtered_df_depth2[which(filtered_df_depth2$sample == i), "SBS_DNMs_Count"] <- 
    nrow(filtered_df_depth2[which(filtered_df_depth2$sample == i), ])
}

# Calculate VAF and coverage metrics
for (sample in unique(filtered_df_depth2$sample)) {
  filtered_df_depth2[which(filtered_df_depth2$sample == sample), "sample.median.vaf"] <- 
    median(filtered_df_depth2[which(filtered_df_depth2$sample == sample), "t_alt_count"] /
             filtered_df_depth2[which(filtered_df_depth2$sample == sample), "t_depth"])
  filtered_df_depth2[which(filtered_df_depth2$sample == sample), "sample.mean.vaf"] <- 
    mean(filtered_df_depth2[which(filtered_df_depth2$sample == sample), "t_alt_count"] /
           filtered_df_depth2[which(filtered_df_depth2$sample == sample), "t_depth"])
  filtered_df_depth2[which(filtered_df_depth2$sample == sample), "Mean_Cov_variants"] <- 
    median(filtered_df_depth2[which(filtered_df_depth2$sample == sample), "t_depth"])
}

# Subset per-sample summary
filtered_df_depth2_per_sample <- unique(filtered_df_depth2[, c("sample", "FamilyID", "FatherID", "MotherID", 
                                                               "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE", 
                                                               "SBS_DNMs_Count", "sample.median.vaf", 
                                                               "sample.mean.vaf", "Mean_Cov_variants")])

# Write filtered table
write.table(filtered_df_depth2_per_sample, 
            file="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_per_offspring_depth2.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

# Define colors based on FamilyID
nb.cols <- length(unique(filtered_df_depth2_per_sample$FamilyID))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
names(mycolors) <- unique(filtered_df_depth2_per_sample$FamilyID)

# Define colors for samples (for VAF plot)
nb.cols_samples <- length(unique(filtered_df_depth2$sample))
sample_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols_samples)
names(sample_colors) <- unique(filtered_df_depth2$sample)

# Plot 1: VAF Distribution (as histogram, color by sample, facet by FamilyID)
vaf_plot_depth2 <- ggplot(filtered_df_depth2, aes(x=t_alt_count/t_depth)) + 
  geom_histogram(aes(fill=sample), alpha=0.4, bins=30, position="identity") + 
  facet_wrap(vars(FamilyID), ncol=2, scales = "free_x") +  # Force each facet to have its own axis
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=14, margin=margin(b=10, t=10)),
        axis.text.y=element_text(size=11), 
        axis.text.x=element_text(angle=0, size=11, vjust=0.5),  # Simplified axis text positioning
        axis.title.x=element_text(size=12, margin=margin(t=30)),  # Increased top margin
        axis.title.y=element_text(size=14),
        plot.margin=unit(c(1, 1, 5, 1), "lines"),  # Significantly increased bottom margin
        plot.subtitle=element_text(hjust=0.5, size=14, margin=margin(t=10, b=10)),
        panel.spacing=unit(3, "lines"),  # Increased panel spacing
        # Remove custom axis line settings temporarily
        # axis.line.x = element_line(color="black", linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray95"),
        plot.title=element_text(hjust=0.5),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) + 
  labs(x="Variant Allele Frequency (VAF)", y="Variant Count", subtitle="Variant Allele Frequency (VAF) Distribution at De Novo Sites") +
  scale_fill_manual(name="Sample", values=sample_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), position = "bottom") +  # Force x-axis position
  coord_cartesian(clip = "off", expand = TRUE)  # Ensure no clipping

# Plot 2: DNM Count by Sample
dnm_by_sample_depth2 <- ggplot(filtered_df_depth2_per_sample, 
                               aes(x=sample, y=SBS_DNMs_Count, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  ylim(0, max(filtered_df_depth2_per_sample$SBS_DNMs_Count) + 50) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Sample") + 
  ylab("De novo variants") + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text.x=element_text(angle=45, hjust=1, size=11),
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Plot 3: DNM Count vs. Mean Coverage
dnm_vs_coverage_depth2 <- ggplot(filtered_df_depth2_per_sample, 
                                 aes(x=MEAN_COVERAGE, y=SBS_DNMs_Count, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, max(filtered_df_depth2_per_sample$SBS_DNMs_Count) + 50) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Mean Coverage") + 
  ylab("De novo variants") + 
  geom_smooth(data=subset(filtered_df_depth2_per_sample, 
                          FamilyID %in% names(which(table(filtered_df_depth2_per_sample$FamilyID) >= 2))),
              aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.line.x=element_line(), 
        axis.line.y=element_line(),
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Plot 4: Distribution of Coverage for De Novo Variants (as histogram)
coverage_dist_depth2 <- ggplot(filtered_df_depth2, aes(x=t_depth)) + 
  geom_histogram(aes(fill=FamilyID), alpha=0.4, bins=30, position="identity") + 
  facet_wrap(vars(sample), ncol=2) + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=12, margin=margin(b=10, t=10)),
        axis.text=element_text(size=10), 
        axis.text.x=element_text(margin=margin(t=10)),
        axis.title.x=element_text(size=12, margin=margin(t=20)), 
        axis.title.y=element_text(size=12),
        plot.margin=unit(c(1, 1, 4, 1), "lines"),
        plot.title=element_text(hjust=0.5),
        panel.spacing=unit(2, "lines"),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray95")) +  
  xlab("Coverage (Depth) at De Novo Sites") + 
  ylab("Variant Count") + 
  scale_fill_manual(name="FamilyID", values=mycolors) + 
  xlim(0, 100) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  facet_wrap(vars(sample), ncol=2, scales = "free_x")  # Add free x scales

# Plot 5: Average Coverage vs. Average VAF
cov_vaf_depth2 <- ggplot(filtered_df_depth2_per_sample, 
                         aes(x=Mean_Cov_variants, y=sample.mean.vaf, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, 1) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Average Coverage at De Novo Sites") + 
  ylab("Average VAF") + 
  geom_smooth(data=subset(filtered_df_depth2_per_sample, 
                          FamilyID %in% names(which(table(filtered_df_depth2_per_sample$FamilyID) >= 2))),
              aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text=element_text(size=10), 
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Save all plots (increased width for VAF plot to accommodate legend)
ggsave(vaf_plot_depth2, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_VAF_distribution_depth2.png", 
       dpi=300, units="in", width=10, height=10)  # Increased width
ggsave(dnm_by_sample_depth2, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_by_Sample_depth2.png", 
       dpi=300, units="in", width=7, height=5)
ggsave(dnm_vs_coverage_depth2, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_vs_Coverage_depth2.png", 
       dpi=300, units="in", width=7, height=5)
ggsave(coverage_dist_depth2, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Coverage_Distribution_DeNovo_depth2.png", 
       dpi=300, units="in", width=7, height=10)
ggsave(cov_vaf_depth2, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Avg_Coverage_vs_Avg_VAF_depth2.png", 
       dpi=300, units="in", width=7, height=5)

#filteredfiledirectlyfilteredinvcf:
#filtered file:
library(ggplot2)
library(RColorBrewer)

# Load filtered data
filtered_df_depth3 <- read.table("D:/PhD/rp793_Cagan/rotationresults/all_denovo_snvs_filtered.txt", 
                                 header=FALSE, sep="\t", 
                                 col.names=c("sample", "chrom", "pos", "ref", "alt", "t_alt_count", "t_depth"))

# Load metadata
metadata <- read.delim("D:/PhD/rp793_Cagan/rotationresults/metadata.txt", 
                       header=TRUE, sep="\t")

# Merge with metadata, keeping only necessary columns
filtered_df_depth3 <- merge.data.frame(filtered_df_depth3, 
                                       metadata[, c("Sample", "FamilyID", "FatherID", "MotherID", "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE")], 
                                       by.x="sample", by.y="Sample")

# Calculate DNM count per sample
for (i in unique(filtered_df_depth3$sample)) {
  filtered_df_depth3[which(filtered_df_depth3$sample == i), "SBS_DNMs_Count"] <- 
    nrow(filtered_df_depth2[which(filtered_df_depth2$sample == i), ])
}

# Calculate VAF and coverage metrics
for (sample in unique(filtered_df_depth3$sample)) {
  filtered_df_depth3[which(filtered_df_depth3$sample == sample), "sample.median.vaf"] <- 
    median(filtered_df_depth3[which(filtered_df_depth3$sample == sample), "t_alt_count"] /
             filtered_df_depth3[which(filtered_df_depth3$sample == sample), "t_depth"])
  filtered_df_depth3[which(filtered_df_depth3$sample == sample), "sample.mean.vaf"] <- 
    mean(filtered_df_depth3[which(filtered_df_depth3$sample == sample), "t_alt_count"] /
           filtered_df_depth3[which(filtered_df_depth3$sample == sample), "t_depth"])
  filtered_df_depth3[which(filtered_df_depth3$sample == sample), "Mean_Cov_variants"] <- 
    median(filtered_df_depth3[which(filtered_df_depth3$sample == sample), "t_depth"])
}

# Subset per-sample summary
filtered_df_depth3_per_sample <- unique(filtered_df_depth3[, c("sample", "FamilyID", "FatherID", "MotherID", 
                                                               "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE", 
                                                               "SBS_DNMs_Count", "sample.median.vaf", 
                                                               "sample.mean.vaf", "Mean_Cov_variants")])

# Write filtered table
write.table(filtered_df_depth3_per_sample, 
            file="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_per_offspring_depth3.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

# Define colors based on FamilyID
nb.cols <- length(unique(filtered_df_depth2_per_sample$FamilyID))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
names(mycolors) <- unique(filtered_df_depth2_per_sample$FamilyID)

# Define colors for samples (for VAF plot)
nb.cols_samples <- length(unique(filtered_df_depth2$sample))
sample_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols_samples)
names(sample_colors) <- unique(filtered_df_depth2$sample)

# Plot 1: VAF Distribution (as histogram, color by sample, facet by FamilyID)
vaf_plot_depth3 <- ggplot(filtered_df_depth3, aes(x=t_alt_count/t_depth)) + 
  geom_histogram(aes(fill=sample), alpha=0.4, bins=30, position="identity") + 
  facet_wrap(vars(FamilyID), ncol=2, scales = "free_x") +  # Force each facet to have its own axis
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=14, margin=margin(b=10, t=10)),
        axis.text.y=element_text(size=11), 
        axis.text.x=element_text(angle=0, size=11, vjust=0.5),  # Simplified axis text positioning
        axis.title.x=element_text(size=12, margin=margin(t=30)),  # Increased top margin
        axis.title.y=element_text(size=14),
        plot.margin=unit(c(1, 1, 5, 1), "lines"),  # Significantly increased bottom margin
        plot.subtitle=element_text(hjust=0.5, size=14, margin=margin(t=10, b=10)),
        panel.spacing=unit(3, "lines"),  # Increased panel spacing
        # Remove custom axis line settings temporarily
        # axis.line.x = element_line(color="black", linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray95"),
        plot.title=element_text(hjust=0.5),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) + 
  labs(x="Variant Allele Frequency (VAF)", y="Variant Count", subtitle="Variant Allele Frequency (VAF) Distribution at De Novo Sites") +
  scale_fill_manual(name="Sample", values=sample_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), position = "bottom") +  # Force x-axis position
  coord_cartesian(clip = "off", expand = TRUE)  # Ensure no clipping

# Plot 2: DNM Count by Sample
dnm_by_sample_depth3 <- ggplot(filtered_df_depth3_per_sample, 
                               aes(x=sample, y=SBS_DNMs_Count, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  ylim(0, max(filtered_df_depth2_per_sample$SBS_DNMs_Count) + 50) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Sample") + 
  ylab("De novo variants") + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text.x=element_text(angle=45, hjust=1, size=11),
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Plot 3: DNM Count vs. Mean Coverage
dnm_vs_coverage_depth3 <- ggplot(filtered_df_depth3_per_sample, 
                                 aes(x=MEAN_COVERAGE, y=SBS_DNMs_Count, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, max(filtered_df_depth2_per_sample$SBS_DNMs_Count) + 50) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Mean Coverage") + 
  ylab("De novo variants") + 
  geom_smooth(data=subset(filtered_df_depth2_per_sample, 
                          FamilyID %in% names(which(table(filtered_df_depth2_per_sample$FamilyID) >= 2))),
              aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.line.x=element_line(), 
        axis.line.y=element_line(),
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Plot 4: Distribution of Coverage for De Novo Variants (as histogram)
coverage_dist_depth3 <- ggplot(filtered_df_depth3, aes(x=t_depth)) + 
  geom_histogram(aes(fill=FamilyID), alpha=0.4, bins=30, position="identity") + 
  facet_wrap(vars(sample), ncol=2) + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=12, margin=margin(b=10, t=10)),
        axis.text=element_text(size=10), 
        axis.text.x=element_text(margin=margin(t=10)),
        axis.title.x=element_text(size=12, margin=margin(t=20)), 
        axis.title.y=element_text(size=12),
        plot.margin=unit(c(1, 1, 4, 1), "lines"),
        plot.title=element_text(hjust=0.5),
        panel.spacing=unit(2, "lines"),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray95")) +  
  xlab("Coverage (Depth) at De Novo Sites") + 
  ylab("Variant Count") + 
  scale_fill_manual(name="FamilyID", values=mycolors) + 
  xlim(0, 100) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  facet_wrap(vars(sample), ncol=2, scales = "free_x")  # Add free x scales

# Plot 5: Average Coverage vs. Average VAF
cov_vaf_depth3 <- ggplot(filtered_df_depth3_per_sample, 
                         aes(x=Mean_Cov_variants, y=sample.mean.vaf, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, 1) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Average Coverage at De Novo Sites") + 
  ylab("Average VAF") + 
  geom_smooth(data=subset(filtered_df_depth3_per_sample, 
                          FamilyID %in% names(which(table(filtered_df_depth3_per_sample$FamilyID) >= 2))),
              aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text=element_text(size=10), 
        axis.title.x=element_text(size=12, margin=margin(t=15)),
        plot.margin=unit(c(1, 1, 3, 1), "lines"))

# Save all plots (increased width for VAF plot to accommodate legend)
ggsave(vaf_plot_depth3, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_VAF_distribution_depth3.png", 
       dpi=300, units="in", width=10, height=10)  # Increased width
ggsave(dnm_by_sample_depth3, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_by_Sample_depth3.png", 
       dpi=300, units="in", width=7, height=5)
ggsave(dnm_vs_coverage_depth3, 
       filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_vs_Coverage_depth3.png", 
       dpi=300, units="in", width=7, height=5)
ggsave(coverage_dist_depth3, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Coverage_Distribution_DeNovo_depth3.png", 
       dpi=300, units="in", width=7, height=10)
ggsave(cov_vaf_depth3, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Avg_Coverage_vs_Avg_VAF_depth3.png", 
       dpi=300, units="in", width=7, height=5)
