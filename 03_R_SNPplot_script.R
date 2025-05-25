# R_SNPplot for RStudio and GitHub

# Load required libraries
library(ggplot2)
library(openxlsx)
library(data.table)

# Prompt for input file
# input_file <- file.choose()  # Optional interactive file selector

# Optional input file
input_filename <- "11v95 snps.txt"
snp_data <- read.table(input_filename, header = TRUE, sep = "\t")

# === Extract ID and set up folder ===
seq_id <- sub("[_ ].*", "", input_filename)  # e.g., "11v95"
today_nice <- format(Sys.Date(), "%d%b%Y")   # e.g., "13May2025"
prefix <- paste0(seq_id, "_SNPplots_", today_nice)  # e.g., "11v95_SNPplots_13May2025"

# === Create output folder ===
output_folder <- file.path(getwd(), prefix)
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# === Clean the data ===
snp_data_clean <- snp_data[!grepl("N", snp_data$`SNP.pattern`), ]
snp_data_clean$SNP_position <- as.numeric(snp_data_clean$sequence_1_PosInContg)

# === Save cleaned files ===
write.table(snp_data_clean, paste0(prefix, "_cleaned_snp_data.txt"), sep = "\t", row.names = FALSE)
write.csv(snp_data_clean, paste0(prefix, "_cleaned_snp_data.csv"), row.names = FALSE)

# === SNP binning ===
bin_width <- 1000
bins <- seq(min(snp_data_clean$SNP_position), max(snp_data_clean$SNP_position), by = bin_width)
binned_snps <- data.frame(Bin = bins[-length(bins)], SNPs_kb = 0)

for (i in seq_along(bins)[-length(bins)]) {
  snps_in_bin <- sum(snp_data_clean$SNP_position >= bins[i] &
                       snp_data_clean$SNP_position < bins[i + 1])
  binned_snps$SNPs_kb[i] <- snps_in_bin / (bin_width / 1000)
}

# === Plotting ===
avg_snp_density <- nrow(snp_data_clean) / 
  (max(snp_data_clean$SNP_position) - min(snp_data_clean$SNP_position)) * 1000
avg_snp_label <- sprintf("%.2f SNPs/kb (%d SNPs)", avg_snp_density, nrow(snp_data_clean))
y_max <- max(binned_snps$SNPs_kb)

snp_plot <- ggplot(snp_data_clean, aes(x = SNP_position)) +
  geom_histogram(aes(y = after_stat(count) / (bin_width / 1000)),
                 binwidth = bin_width, fill = "black", alpha = 0.6, color = "black") +
  geom_hline(yintercept = avg_snp_density, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = mean(snp_data_clean$SNP_position), y = y_max * 0.95,
           label = avg_snp_label, color = "red", size = 4, fontface = "bold") +
  labs(title = paste("SNP Density:", seq_id),
       x = "SNP Positions (bp)", y = "SNPs per kb") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  )

# Show the plot
print(snp_plot)

# === Save plot as PNG and PDF ===
ggsave(paste0(prefix, "_snp_density_plot[scale1.5].png"), plot = snp_plot,
       width = 10, height = 4, dpi = 300, units = "in", scale = 1.5)
ggsave(paste0(prefix, "_snp_density_plot[scale1.5].pdf"), plot = snp_plot,
       width = 10, height = 4, units = "in", scale = 1.5)

# === Save SNP bins ===
write.table(binned_snps, paste0(prefix, "_snp_bins_output.txt"), sep = "\t", row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "SNPs per kb")
writeData(wb, "SNPs per kb", binned_snps)
saveWorkbook(wb, paste0(prefix, "_snp_bins_output.xlsx"), overwrite = TRUE)

# Done
message("✅ Analysis complete. Results saved to: ", output_folder)
