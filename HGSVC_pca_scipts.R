# Load required libraries
library(dbscan)
library(fpc)
library(factoextra)
library(SNPRelate)
library(data.table)
library(raster)

# Read in the 1KG info file
hgsvc_samples <- read.delim("1KG_3202_samples.ped", header = TRUE)


# Uncomment the following if the 2504-sample file is needed
# hgsvc_samples_2504 <- read.table("all_2504_1KG_samples.txt", header = FALSE)
# colnames(hgsvc_samples_2504) <- "SampleID"

##### IGNORE THIS BLOCK -- RELEVANT TO HPRC AND HGSVC SAMPLES ONLY
# Read subsets of samples of interest
year_1_samples <- read.table("YEAR_1_LongRead_samples_v2.txt", header = TRUE)
year_2_samples <- read.delim("YEAR_2_LongRead_samples_v2.txt", header = TRUE)
colnames(year_1_samples) <- "Sample"
colnames(year_2_samples) <- "Sample"
hgsvc_lr_YR1and2_samples <- unique(rbind(year_1_samples, year_2_samples))

hprc_samples <- read.delim("HPRC_combined_120samples.txt", header = TRUE)
# (Note: adjust the column name if needed; here we set the second column to 'Sample')
colnames(hprc_samples)[2] <- "Sample"

# Merge HGSVC samples (using the 2504 file if applicable)
merged_1 <- merge(hgsvc_samples, hgsvc_samples_2504, by = "SampleID", all = TRUE)
write.table(merged_1, "merged_1.tab", quote = FALSE, sep = "\t", row.names = FALSE)

# Create a combined list of HGSVC and HPRC sample IDs
hgsvc_hprs_sampleIDs <- unique(sort(c(year_1_samples$Sample, year_2_samples$Sample, hprc_samples$Sample)))



#### PCA ANALYSIS ####

vcf.fn <- "pangenie_merged_bi_nosnvs.vcf.gz"

# Convert VCF to GDS format and open the file
snpgdsVCF2GDS(vcf.fn, "hgsvc_sv.gds", method = "copy.num.of.ref", verbose = TRUE)
genofile <- openfn.gds("hgsvc_sv.gds")

# Run PCA using the randomized algorithm
pca_res <- snpgdsPCA(genofile, num.thread = 4, algorithm = "randomized")
eigenvect <- pca_res$eigenvect

# Build a data frame with PCA results and annotate with 1KG population info
out_mx <- data.frame(SampleID = pca_res$sample.id, eigenvect)
n_pcs <- ncol(out_mx) - 1
colnames(out_mx)[-1] <- paste0("PC", 1:n_pcs)
merged_out <- merge(out_mx, hgsvc_samples, by = "SampleID")
write.table(merged_out, "HGSVC_3202_Autosomes_PCA.tab", quote = FALSE, sep = "\t", row.names = FALSE)

# If rerunning, load the saved PCA matrix (adjust file path as needed)
merged_out <- read.delim("HGSVC_3202_Autosomes_PCA.tab", header = TRUE)
# Optionally, to use a different PCA result:
merged_out <- read.delim("SNV_genotypes/HGSVC_2504_chr22GATK_PCA.tab", header = TRUE)

# Coerce the first 7 PCs to numeric (you can adjust the range if more PCs are present)
for(i in 1:7){
  colname <- paste0("PC", i)
  merged_out[[colname]] <- as.numeric(as.character(merged_out[[colname]]))
}

# Set default color and override for each superpopulation
merged_out$Superpopulation_colors <- "black"
merged_out$Superpopulation_colors[merged_out$Superpopulation == "EUR"] <- "#e41a1c"
merged_out$Superpopulation_colors[merged_out$Superpopulation == "EAS"] <- "#377eb8"
merged_out$Superpopulation_colors[merged_out$Superpopulation == "SAS"] <- "#4daf4a"
merged_out$Superpopulation_colors[merged_out$Superpopulation == "AFR"] <- "#984ea3"
merged_out$Superpopulation_colors[merged_out$Superpopulation == "AMR"] <- "#ff7f00"

# Basic PCA plots for Pangenie genotypes
par(mfrow = c(1,1))
plot(merged_out$PC1, merged_out$PC2, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC1", ylab = "PC2")
plot(merged_out$PC2, merged_out$PC3, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC2", ylab = "PC3", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065))

# PCA panels for GATK genotypes
par(mfrow = c(2,2))
plot(merged_out$PC1, merged_out$PC2, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC1 (1.29%)", ylab = "PC2 (0.42%)", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.05))
plot(merged_out$PC2, merged_out$PC3, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC2 (0.42%)", ylab = "PC3 (0.21%)", xlim = c(-0.05, 0.05), ylim = c(-0.05, 0.065))
plot(merged_out$PC3, merged_out$PC4, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC3 (0.21%)", ylab = "PC4 (0.16%)", xlim = c(-0.05, 0.065), ylim = c(-0.04, 0.12))
plot(merged_out$PC4, merged_out$PC5, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC4 (0.16%)", ylab = "PC5 (0.14%)", xlim = c(-0.04, 0.12), ylim = c(-0.25, 0.052))
# (Additional plots for PC5 vs PC6 and PC6 vs PC7 are generated below without mfrow reset)

plot(merged_out$PC5, merged_out$PC6, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC5", ylab = "PC6")
plot(merged_out$PC6, merged_out$PC7, pch = 19, col = merged_out$Superpopulation_colors,
     xlab = "PC6", ylab = "PC7")

# Label specific subsets: YEAR 1 and YEAR 2 samples
merged_out_yr1 <- merged_out[merged_out$SampleID %in% year_1_samples$YR_1_SAMPLE, ]
merged_out_yr2 <- merged_out[merged_out$SampleID %in% year_2_samples$Sample, ]
merged_out_yr1_and_2 <- rbind(merged_out_yr1, merged_out_yr2)

# HPRC subset (note: adjust column name if necessary)
merged_out_hprc <- merged_out[merged_out$SampleID %in% hprc_samples$ChildID, ]
combined_hprc_hgsvc <- rbind(merged_out_hprc, merged_out_yr1_and_2)

# Overlay HPRC samples as filled diamonds
par(new = TRUE)
plot(merged_out_hprc$PC1, merged_out_hprc$PC2, pch = 23, col = merged_out_hprc$Superpopulation_colors,
     xlab = "", ylab = "", bg = "darkgrey", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.03), cex = 1.2)
par(new = TRUE)
plot(merged_out_hprc$PC2, merged_out_hprc$PC3, pch = 23, col = merged_out_hprc$Superpopulation_colors,
     xlab = "", ylab = "", bg = "grey", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065), cex = 1.2)

# Overlay YEAR 1 and 2 Long Read samples as triangles filled in black
par(new = TRUE)
plot(merged_out_yr1_and_2$PC1, merged_out_yr1_and_2$PC2, pch = 25, col = merged_out_yr1_and_2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.03), cex = 1.2)
par(new = TRUE)
plot(merged_out_yr1_and_2$PC2, merged_out_yr1_and_2$PC3, pch = 25, col = merged_out_yr1_and_2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065), cex = 1.2)



# Overlay centroid approach samples (MIN as circles, MAX as squares)
merged_out_centroidA <- merged_out[merged_out$SampleID %in% centroid_merged_optionA$CENTROID_SAMPLES_MIN, ]
merged_out_centroidB <- merged_out[merged_out$SampleID %in% centroid_merged_optionB$CENTROID_SAMPLES_MAX, ]

par(new = TRUE)
plot(merged_out_centroidA$PC1, merged_out_centroidA$PC2, pch = 21, col = merged_out_centroidA$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.03), cex = 1.2)
par(new = TRUE)
plot(merged_out_centroidB$PC1, merged_out_centroidB$PC2, pch = 22, col = merged_out_centroidB$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.03), cex = 1.2)
par(new = TRUE)
plot(merged_out_centroidA$PC2, merged_out_centroidA$PC3, pch = 21, col = merged_out_centroidA$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065), cex = 1.2)
par(new = TRUE)
plot(merged_out_centroidB$PC2, merged_out_centroidB$PC3, pch = 22, col = merged_out_centroidB$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065), cex = 1.2)

# Additional overlay for combined HPRC and HGSVC projections (using GATK genotypes)
par(new = TRUE)
plot(merged_out_yr1_and_2$PC1, merged_out_yr1_and_2$PC2, pch = 25, col = merged_out_yr1_and_2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.05))
par(new = TRUE)
plot(merged_out_yr1_and_2$PC2, merged_out_yr1_and_2$PC3, pch = 25, col = merged_out_yr1_and_2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065))
par(new = TRUE)
plot(merged_out_yr1_and_2$PC3, merged_out_yr1_and_2$PC4, pch = 25, col = merged_out_yr1_and_2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.065), ylim = c(-0.04, 0.12))
par(new = TRUE)
plot(merged_out_yr1_and_2$PC4, merged_out_yr1_and_2$PC5, pch = 25, col = merged_out_yr1_and_2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.04, 0.12), ylim = c(-0.25, 0.052))


# Individual overlays for YEAR 1 and YEAR 2 samples (if needed)
par(new = TRUE)
plot(merged_out_yr1$PC1, merged_out_yr1$PC2, pch = 25, col = merged_out_yr1$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black")
par(new = TRUE)
plot(merged_out_yr1$PC2, merged_out_yr1$PC3, pch = 25, col = merged_out_yr1$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065))
par(new = TRUE)
plot(merged_out_yr2$PC1, merged_out_yr2$PC2, pch = 25, col = merged_out_yr2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.02, 0.05), ylim = c(-0.05, 0.03))
par(new = TRUE)
plot(merged_out_yr2$PC2, merged_out_yr2$PC3, pch = 25, col = merged_out_yr2$Superpopulation_colors,
     xlab = "", ylab = "", bg = "black", xlim = c(-0.05, 0.03), ylim = c(-0.05, 0.065))



#### CLUSTERING APPROACHES ####

# Example: kmeans on EUR samples (using PC2 and PC3)
merged_kmeans <- kmeans(merged_out[merged_out$Superpopulation == "EUR", 2:3], centers = 1)

# Determine the optimal number of clusters for AMR samples using gap statistic
fviz_nbclust(merged_out[merged_out$Superpopulation == "AMR", 2:3], kmeans, method = "gap_stat")

# Perform kmeans clustering with user-defined centers per superpopulation
km.opt.eur <- kmeans(merged_out[merged_out$Superpopulation == "EUR", 2:3], centers = 3, nstart = 25)
km.opt.afr <- kmeans(merged_out[merged_out$Superpopulation == "AFR", 2:3], centers = 7, nstart = 25)
km.opt.eas <- kmeans(merged_out[merged_out$Superpopulation == "EAS", 2:3], centers = 1, nstart = 25)
km.opt.sas <- kmeans(merged_out[merged_out$Superpopulation == "SAS", 2:3], centers = 4, nstart = 25)
km.opt.amr <- kmeans(merged_out[merged_out$Superpopulation == "AMR", 2:3], centers = 6, nstart = 25)

# Visualize one of the clusters (e.g., AMR)
fviz_cluster(km.opt.amr,
             data = merged_out[merged_out$Superpopulation == "AMR", 2:3],
             ellipse.type = "convex", palette = "jco", ggtheme = theme_minimal())



#### DISTANCE FROM CENTROIDS ####

# Prepare a data frame to store distances; adjust columns (1,18:24) as needed
merged_out_distances_from_centroids <- merged_out[, c(1, 18:24)]
merged_out_distances_from_centroids$DIST_from_CENTROID <- rep(99999, nrow(merged_out))
merged_out_distances_from_centroids$CENTROID_ID <- rep(0, nrow(merged_out))

# List of kmeans objects to iterate over
kmeans_array <- c("km.opt.afr", "km.opt.eur", "km.opt.eas", "km.opt.sas", "km.opt.amr")

for(i in seq_along(kmeans_array)){
  current_km <- eval(parse(text = kmeans_array[i]))
  current_clusters <- current_km$cluster
  
  for(j in sort(unique(current_clusters))){
    slice_idx <- as.numeric(names(which(current_clusters == j)))
    center_coords <- current_km$centers[j, ]
    # Create a matrix where each row is the center coordinates
    center_matrix <- matrix(rep(center_coords, length(slice_idx)), ncol = 2, byrow = TRUE)
    # Calculate Euclidean distances from each sample to its cluster center
    distances <- pointDistance(merged_out[slice_idx, 2:3], center_matrix, lonlat = FALSE)
    
    merged_out_distances_from_centroids[slice_idx, "DIST_from_CENTROID"] <- distances
    merged_out_distances_from_centroids[slice_idx, "CENTROID_ID"] <- j
  }
}

write.csv(merged_out_distances_from_centroids,
          "merged_out_distances_from_centroids_21.csv",
          row.names = FALSE, quote = FALSE)

# Reload and filter out the selected samples
merged_out_distances_from_centroids <- read.csv("merged_out_distances_from_centroids_21.csv", header = TRUE)
merged_out_distances_from_centroids <- merged_out_distances_from_centroids[
  !merged_out_distances_from_centroids$SampleID %in% hgsvc_hprs_sampleIDs, ]

# Identify extreme (min and max distance) samples for each superpopulation and cluster
centroid_merged_optionA <- c(NA)
centroid_merged_optionB <- c(NA)

for (superpop in c("AFR", "AMR", "EUR", "EAS", "SAS")) {
  slice <- merged_out_distances_from_centroids[merged_out_distances_from_centroids$Superpopulation == superpop, ]
  
  for (cluster in sort(unique(slice$CENTROID_ID))) {
    slice_cluster <- slice[slice$CENTROID_ID == cluster, ]
    sample_min <- slice_cluster[which.min(slice_cluster$DIST_from_CENTROID), "SampleID"]
    sample_max <- slice_cluster[which.max(slice_cluster$DIST_from_CENTROID), "SampleID"]
    centroid_merged_optionA <- c(centroid_merged_optionA, sample_min)
    centroid_merged_optionB <- c(centroid_merged_optionB, sample_max)
  }
}

centroid_merged_optionA <- data.frame(CENTROID_SAMPLES_MIN = na.omit(as.character(centroid_merged_optionA)))
centroid_merged_optionB <- data.frame(CENTROID_SAMPLES_MAX = na.omit(as.character(centroid_merged_optionB)))
