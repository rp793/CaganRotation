library(ggplot2)

# Load de novo SNVs
all_Samples_snv_dnms <- read.table("D:/PhD/rp793_Cagan/rotationresults/all_denovo_snvs.txt", 
                                   header=FALSE, sep="\t", 
                                   col.names=c("sample", "chrom", "pos", "ref", "alt", "t_alt_count", "t_depth"))

# Load metadata
metadata <- read.delim("D:/PhD/rp793_Cagan/rotationresults/metadata.txt", 
                       header=TRUE, sep="\t")

# Merge with metadata
all_Samples_snv_dnms <- merge.data.frame(all_Samples_snv_dnms, metadata, 
                                         by.x="sample", by.y="Sample")

# Calculate DNM count per sample
for (i in unique(all_Samples_snv_dnms$sample)) {
  count <- nrow(all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == i), ])
  all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == i), "SBS_DNMs_Count"] <- count
  print(paste("Sample:", i, "DNM Count:", count))
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

# Subset for per-sample summary
all_Samples_snv_dnms_per_sample <- unique(all_Samples_snv_dnms[, c("sample", "FamilyID", "FatherID", "MotherID", 
                                                                   "Sex", "MEAN_COVERAGE", "MEDIAN_COVERAGE", 
                                                                   "SBS_DNMs_Count", "sample.median.vaf", 
                                                                   "sample.mean.vaf", "Mean_Cov_variants")])

# Write table
write.table(all_Samples_snv_dnms_per_sample, 
            file="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_per_offspring.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

# VAF distribution plot
VAF_plot <- ggplot(all_Samples_snv_dnms, aes(x = t_alt_count/t_depth)) + 
  geom_density(aes(group=sample), alpha=0.4, fill="grey") + 
  facet_wrap(vars(FamilyID), dir="v", ncol=2) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.x = element_text(hjust=0.5, size=14), 
        axis.text.y = element_text(size=11), axis.text.x = element_text(angle=0, size=11), 
        axis.title = element_text(size=14)) + 
  xlab("VAF") + 
  xlim(0,1)

ggsave(VAF_plot, filename="D:/PhD/rp793_Cagan/rotationresults/SNV_VAF_distribution.png", 
       dpi=300, units="in", width=7, height=10)

library(ggplot2)
library(RColorBrewer)  # For color palette
# library(ggpubr)      # Uncomment if you want theme_pubclean(), otherwise skip

# Assuming this is already run up to write.table()
# Load your data if starting fresh
all_Samples_snv_dnms_per_sample <- read.delim("D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_per_offspring.txt", 
                                              header=TRUE, sep="\t")

# Set y-axis max for plots
y_max <- max(all_Samples_snv_dnms_per_sample$SBS_DNMs_Count) + 50

# Since no age data, skip age sorting and labeling
# Instead, use sample directly (no Tumor_Sample_Barcode_label needed yet)
all_Samples_snv_dnms_per_sample$sample <- factor(all_Samples_snv_dnms_per_sample$sample)

# Define colors based on FamilyID
nb.cols <- length(unique(all_Samples_snv_dnms_per_sample$FamilyID))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
names(mycolors) <- unique(all_Samples_snv_dnms_per_sample$FamilyID)
unname(mycolors)

# Plot 1: DNM Count by Sample, colored by FamilyID (replacing age plot)
p <- ggplot() + 
  geom_point(data=all_Samples_snv_dnms_per_sample, 
             aes(x=sample, y=SBS_DNMs_Count, fill=as.character(FamilyID)), 
             size=3.5, pch=21, stroke=0.25) + 
  ylim(0, y_max) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Sample") + 
  ylab("De novo variants") + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text.x=element_text(angle=45, hjust=1, size=11))  # Rotate x-axis labels

# Plot 2: DNM Count vs MEAN_COVERAGE (alternative to paternal age)
p_BurdenCov_snv <- ggplot() + 
  geom_point(data=all_Samples_snv_dnms_per_sample, 
             aes(x=MEAN_COVERAGE, y=SBS_DNMs_Count, fill=as.character(FamilyID)), 
             size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, y_max) + 
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Mean Coverage") + 
  ylab("De novo variants") + 
  geom_smooth(data=all_Samples_snv_dnms_per_sample, 
              aes(x=MEAN_COVERAGE, y=SBS_DNMs_Count, 
                  group=as.character(FamilyID), 
                  color=as.character(FamilyID)), 
              method='lm', se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() +  # Use theme_bw() instead of theme_pubclean() if ggpubr isn’t installed
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.line.x=element_line(), 
        axis.line.y=element_line())

# Save Plot 1
ggsave(p, filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_by_Sample.png", 
       dpi=300, units="in", width=7, height=5)

# Save Plot 2
pdf(file="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_vs_Coverage.pdf", width=7, height=5)
print(p_BurdenCov_snv)
dev.off()
ggsave(p_BurdenCov_snv, filename="D:/PhD/rp793_Cagan/rotationresults/SNV_DNMs_vs_Coverage.png", 
       dpi=300, units="in", width=7, height=5)

#additional plots:
library(ggplot2)
library(RColorBrewer)  # For color palette

# Load data (assuming already run, or reload if needed)
all_Samples_snv_dnms <- read.table("D:/PhD/rp793_Cagan/rotationresults/all_denovo_snvs.txt", 
                                   header=FALSE, sep="\t", 
                                   col.names=c("sample", "chrom", "pos", "ref", "alt", "t_alt_count", "t_depth"))
metadata <- read.delim("D:/PhD/rp793_Cagan/rotationresults/metadata.txt", header=TRUE, sep="\t")
all_Samples_snv_dnms <- merge.data.frame(all_Samples_snv_dnms, metadata, by.x="sample", by.y="Sample")

# Calculate DNM count (if not already done)
for (i in unique(all_Samples_snv_dnms$sample)) {
  all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == i), "SBS_DNMs_Count"] <- 
    nrow(all_Samples_snv_dnms[which(all_Samples_snv_dnms$sample == i), ])
}

# Calculate VAF and coverage metrics (if not already done)
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

# Define colors based on FamilyID
nb.cols <- length(unique(all_Samples_snv_dnms_per_sample$FamilyID))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
names(mycolors) <- unique(all_Samples_snv_dnms_per_sample$FamilyID)

# Plot 1: Distribution of coverage for all called de novo variants per animal
coverage_dist_plot <- ggplot(all_Samples_snv_dnms, aes(x=t_depth)) + 
  geom_density(aes(group=sample, fill=FamilyID), alpha=0.4) + 
  facet_wrap(vars(sample), ncol=2) + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(hjust=0.5, size=12), 
        axis.text=element_text(size=10), 
        axis.title=element_text(size=12)) + 
  xlab("Coverage (Depth) at De Novo Sites") + 
  ylab("Density") + 
  scale_fill_manual(values=mycolors) +
  xlim(0, 100)  # Adjust range based on your data if needed

# Plot 2: Average coverage vs. average VAF per animal
cov_vaf_plot <- ggplot(all_Samples_snv_dnms_per_sample, 
                       aes(x=Mean_Cov_variants, y=sample.mean.vaf, fill=FamilyID)) + 
  geom_point(size=3.5, pch=21, stroke=0.25) + 
  xlim(0, NA) + 
  ylim(0, 1) +  # VAF ranges 0–1
  scale_fill_manual(name="Family", values=mycolors) + 
  xlab("Average Coverage at De Novo Sites") + 
  ylab("Average VAF") + 
  geom_smooth(aes(group=FamilyID, color=FamilyID), method="lm", se=FALSE) + 
  scale_color_manual(name="Family", values=mycolors) + 
  theme_bw() + 
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10), 
        legend.background=element_blank(), 
        legend.position="top", 
        legend.direction="horizontal",
        axis.text=element_text(size=10), 
        axis.title=element_text(size=12))

# Save Plot 1
ggsave(coverage_dist_plot, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Coverage_Distribution_DeNovo.png", 
       dpi=300, units="in", width=7, height=10)

# Save Plot 2
ggsave(cov_vaf_plot, 
       filename="D:/PhD/rp793_Cagan/rotationresults/Avg_Coverage_vs_Avg_VAF.png", 
       dpi=300, units="in", width=7, height=5)

#edited to improve plotting on original data: 
