library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ape)

beta = read.csv("/home/rodrigo/drought_results/postprocessing/beta-div.tsv", sep="\t")
beta$comparison1 <- as.character(beta$comparison1)
beta$comparison2 <- as.character(beta$comparison2)
bc_mat <- beta[,c(2,3,4)]

# Get unique sample names
samples <- unique(c(bc_mat$comparison1, bc_mat$comparison2))

# Create an empty matrix with NAs
triangular_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))

# Fill in the matrix with the distance values
for (i in 1:nrow(bc_mat)) {
  row_name <- bc_mat$comparison1[i]
  col_name <- bc_mat$comparison2[i]
  value <- bc_mat$braycurtis[i]
  triangular_matrix[row_name, col_name] <- value
}

# Replace missing values (diagonal and lower triangle) with 0
triangular_matrix[is.na(triangular_matrix) | lower.tri(triangular_matrix)] <- 0
#print(triangular_matrix)

pcoa <- as.data.frame(cmdscale(triangular_matrix, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

my_labels = gsub("[_].*", "", rownames(pcoa))
plot(pcoa$PCoA1, pcoa$PCoA2, col=as.factor(my_labels), xlab="PCoA1", ylab="PCoA2")
legend("top", unique(my_labels), pch=1, col=(unique(as.factor(my_labels))))
#text(pcoa$PCoA1, pcoa$PCoA2, samples, pos=3)





beta_data <- read.delim("/home/rodrigo/drought_results/postprocessing/beta-div.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(c(beta_data$comparison1, beta_data$comparison2))
dist_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples),
                      dimnames = list(samples, samples))

for (i in 1:nrow(beta_data)) {
    s1 <- beta_data$comparison1[i]
    s2 <- beta_data$comparison2[i]
    bc <- beta_data$braycurtis[i]
    dist_matrix[s1, s2] <- bc
    dist_matrix[s2, s1] <- bc
}
diag(dist_matrix) <- 0
dist_object <- as.dist(dist_matrix)

metadata_sra <- read.csv("/home/rodrigo/Downloads/SraRunTable2.csv", header = TRUE, stringsAsFactors = FALSE)
metadata_sra$Sample <- trimws(metadata_sra$Run)

metadata_lake <- read.csv("/home/rodrigo/Downloads/amostras_lago.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(metadata_lake)[colnames(metadata_lake) == "Lake"] <- "Lake_of_Origin"

combined_metadata <- merge(metadata_sra, metadata_lake, by = "Sample", all.x = TRUE)

combined_metadata$Host_diet_Clean <- trimws(as.character(combined_metadata$Host_diet))
combined_metadata$Lake_of_Origin_Clean <- trimws(as.character(combined_metadata$Lake_of_Origin))

combined_metadata$Grupo_Final <- paste0(combined_metadata$Host_diet_Clean, "_", combined_metadata$Lake_of_Origin_Clean)

samples_in_dist_matrix <- rownames(as.matrix(dist_object))

filtered_metadata_present_in_dist <- combined_metadata %>%
  filter(Sample %in% samples_in_dist_matrix)

final_merged_data <- filtered_metadata_present_in_dist %>%
  filter(Host_diet_Clean %in% c("benthic", "limnetic")) %>%
  filter(!is.na(Grupo_Final) & Grupo_Final != "")

samples_for_nmds <- final_merged_data %>% pull(Sample)

dist_subset <- as.dist(as.matrix(dist_object)[samples_for_nmds, samples_for_nmds])

desired_group_order <- c(
  "benthic_Paxton", "benthic_Priest", "benthic_Little Quarry",
  "limnetic_Paxton", "limnetic_Priest", "limnetic_Little Quarry"
)
actual_groups_present <- intersect(desired_group_order, unique(final_merged_data$Grupo_Final))
final_merged_data$Grupo_Final <- factor(final_merged_data$Grupo_Final, levels = actual_groups_present)

nmds_result <- metaMDS(dist_subset, autotransform = FALSE, wascores = FALSE, k = 2, trymax = 100)

print(paste("NMDS Stress:", nmds_result$stress))
plot(nmds_result, type = "t")

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)

final_nmds_data_for_plot <- merge(nmds_scores, final_merged_data[, c("Sample", "Grupo_Final", "Host_diet_Clean")], by = "Sample", all.x = TRUE)

final_nmds_data_for_plot <- final_nmds_data_for_plot %>%
    filter(!is.na(Grupo_Final))

hulls <- final_nmds_data_for_plot %>%
    group_by(Grupo_Final) %>%
    slice(chull(NMDS1, NMDS2))

num_unique_groups <- length(unique(final_nmds_data_for_plot$Grupo_Final))

colors_for_groups <- brewer.pal(min(num_unique_groups, 8), "Dark2")
if (num_unique_groups > 8) {
    colors_for_groups <- c(brewer.pal(8, "Dark2"), brewer.pal(min(num_unique_groups - 8, 8), "Set2"))
}

shapes_for_groups <- c(16, 17, 15, 18, 8, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12)[1:num_unique_groups]

ggplot(final_nmds_data_for_plot, aes(x = NMDS1, y = NMDS2, color = Grupo_Final, fill = Grupo_Final)) +
    geom_point(size = 3, aes(shape = Grupo_Final)) +
    geom_polygon(data = hulls, alpha = 0.3, show.legend = FALSE) +
    scale_color_manual(values = colors_for_groups) +
    scale_shape_manual(values = shapes_for_groups) +
    labs(title = "NMDS da Composição do Microbioma por Dieta e Lago (Benthic e Limnetic)",
         x = "NMDS1",
         y = "NMDS2",
         color = "Grupo (Dieta_Lago)",
         fill = "Grupo (Dieta_Lago)") +
    theme_minimal() +
    theme(legend.position = "right")
