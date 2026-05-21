install.packages("dplyr")
library(dplyr)
install.packages("BiocManager")
BiocManager::install("pwalign")
BiocManager::install("Biostrings")
library(Biostrings)
library(pwalign)









set.seed(123)

base_counts <- matrix(
  rpois(5000 * 6, lambda = 80),
  nrow = 5000,
  ncol = 6
)

depth_factors <- c(0.7, 1.0, 1.3, 0.8, 1.2, 1.5)

counts <- sweep(base_counts, 2, depth_factors, "*")

counts <- round(counts)

rownames(counts) <- paste0("Gene_", 1:5000)

colnames(counts) <- c(
  "Control_1",
  "Control_2",
  "Control_3",
  "Treated_1",
  "Treated_2",
  "Treated_3"
)

write.csv(counts, "counts.csv")

head(counts)
dim(counts)




#PART 1#
counts <- read.csv("C:/Users/rahaf/OneDrive/سطح المكتب/BIOOOO/counts (2).csv",row.names = 1)
head(counts)
gene_totals <- rowSums(counts)
head(gene_totals)

colnames(counts) <- c(
  "Control_1",
  "Control_2",
  "Control_3",
  "Treated_1",
  "Treated_2",
  "Treated_3"
)
str(counts)
dim(counts)
summary(counts)

filtered_counts <- counts[gene_totals >= 15, ]
dim(filtered_counts)

set.seed(123)

n_rows <- nrow(filtered_counts)
n_cols <- ncol(filtered_counts)

random_rows <- sample(1:n_rows, 5)
random_cols <- sample(1:n_cols, 5, replace = TRUE)

for(i in 1:5){
  filtered_counts[random_rows[i], random_cols[i]] <- NA
}
na_locations <- which(is.na(filtered_counts), arr.ind = TRUE)

print("Locations of NA values:")
print(na_locations)

for(i in 1:nrow(na_locations)){
  
  row <- na_locations[i, 1]
  col <- na_locations[i, 2]
  
  if(col <= 3){
    
    values <- as.numeric(filtered_counts[row, 1:3])
    
    replacement_value <- sum(values, na.rm = TRUE) / sum(!is.na(values))
    
  } else {
    
    values <- as.numeric(filtered_counts[row, 4:6])
    replacement_value <- sum(values, na.rm = TRUE) / sum(!is.na(values))
    
  }
  
  filtered_counts[row, col] <- replacement_value
}

sum(is.na(filtered_counts))

sample_totals <- colSums(filtered_counts)
print(sample_totals)

cpm_matrix <- (filtered_counts / sample_totals) * 1000000

head(cpm_matrix)

log_counts <- log2(filtered_counts + 1)

par(mfrow = c(1,2))

boxplot(log_counts,
        main = "Log2(counts + 1)",
        col = "lightblue")
boxplot(cpm_matrix,
        main = "CPM Normalized Counts",
        col = "lightgreen")




#PART 2#

wt_seq <- readDNAStringSet("C:/Users/rahaf/OneDrive/سطح المكتب/BIOOOO/TP53_WT (2).fasta")[[1]]
mut_seq <- readDNAStringSet("C:/Users/rahaf/OneDrive/سطح المكتب/BIOOOO/TP53_Variant (2).fasta")[[1]]
wt_seq
mut_seq
alignment <- pwalign::pairwiseAlignment(
  wt_seq,
  mut_seq,
  type = "global"
)
alignment
alignment_score <- score(alignment)
alignment_score



wt_rna <- chartr("T", "U", as.character(wt_seq))
wt_rna <- RNAString(wt_rna)
wt_rna
mut_rna <- chartr("T", "U", as.character(mut_seq))
mut_rna <- RNAString(mut_rna)
mut_rna


variant_protein <- translate(wt_rna)
variant_protein
variant_protein_mut <- translate(mut_rna)
variant_protein



gc_content <- letterFrequency(wt_seq, letters = "GC", as.prob = TRUE) * 100
gc_content

protein_char <- as.character(variant_protein_mut)

star_positions <- gregexpr("\\*", protein_char)[[1]]
protein_length <- nchar(protein_char)

if (star_positions[1] == -1) {
  print("No STOP codon detected\n")
} else if (protein_length %in% star_positions) {
  print("Full Length Protein (Normal STOP codon at the end)\n")
} else {
  print("Truncated Protein Detected (Premature STOP codon)\n")
}




#PART 3#



set.seed(456)
prot_intensities <- matrix(rnorm(4800, mean = 20, sd = 4), ncol = 6)

rownames(prot_intensities) <- paste0("Prot_", 1:800)
colnames(prot_intensities) <- c("C1", "C2", "C3", "T1", "T2", "T3")

prot_intensities[sample(1:4800, 400)] <- NA
 
total_missing <- sum(is.na(prot_intensities))
print(paste("Total missing values:", total_missing))

filter_proteins <- function(x) {
  ctrl_valid <- sum(!is.na(x[1:3])) >= 2
  treat_valid <- sum(!is.na(x[4:6])) >= 2
  return(ctrl_valid | treat_valid)}
keep_idx <- apply(prot_intensities, 1, filter_proteins)
filtered_prot <- prot_intensities[keep_idx, ]

min_val <- min(filtered_prot, na.rm = TRUE)
print(paste("Minimum detected intensity value:", min_val))









#P4
min_value <- min(filtered_prot, na.rm = TRUE)
imputation_value <- min_value / 5

min_value
imputation_value

filtered_prot[is.na(filtered_prot)] <- imputation_value

sum(is.na(filtered_prot))  


groupA_mean <- rowMeans(filtered_prot[, 1:3])
groupB_mean <- rowMeans(filtered_prot[, 4:6])

groupA_mean
groupB_mean


fold_change <- groupB_mean / groupA_mean

fold_change


z_score_prot <- t(apply(filtered_prot, 1, function(x) {
  (x - mean(x)) / sd(x)
}))


protein_df <- data.frame(
  GroupA = groupA_mean,
  GroupB = groupB_mean,
  FoldChange = fold_change
)

protein_df


viz_df <- protein_df


library(ggplot2)


ggplot(viz_df, aes(x = GroupA, y = GroupB, color = FoldChange)) +
  geom_point() +
  labs(
    title = "Protein Abundance: Group A vs Group B",
    x = "Group A Mean Intensity",
    y = "Group B Mean Intensity",
    color = "Fold Change"
  )


upregulated <- protein_df$FoldChange > 1
downregulated <- protein_df$FoldChange < 1

upregulated
downregulated


