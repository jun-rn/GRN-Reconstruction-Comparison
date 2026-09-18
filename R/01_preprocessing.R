#extracting the .tar file
untar("GSE324044_RAW.tar", exdir = "GSE324044_samples/" )

library(tidyverse)

rsem_files <- list.files("GSE324044_samples/",
                         pattern = "\\.genes\\.results\\.gz$",
                         full.names = TRUE)

#check if files are found
print(rsem_files)

#read and merge all samples into one matrix
merged <- rsem_files %>%
  lapply(function(f){
    #extract sample name from filename
    sample_name <- basename(f) %>%
      str_remove("\\.genes\\.results\\.gz$")
    
    #read the file
    read_tsv(f, show_col_types = FALSE) %>%
      select(gene_id, TPM) %>%
      rename(!!sample_name := TPM)
  }) %>%
  reduce(full_join, by = "gene_id")
#set gene_id as row names
merged <- merged %>% column_to_rownames("gene_id")

#Quick Check
dim(merged)  # should print: number of genes x number of samples
head(merged[, 1:5])  # preview first 5 samples

#Filter Low Expression Genes
#Keep gene only if TPM>1 in atleast 20% of samples
min_samples <- ceiling(0.2 * ncol(merged))

filtered <- merged[rowSums(merged > 1) >= min_samples, ]

cat("Genes before filtering:", nrow(merged), "\n")
cat("Genes after filtering:", nrow(filtered), "\n")

#Convert to HUGO Names
library(biomaRt)

ensembl_ids <- rownames(filtered) %>% str_remove("\\..*$")

# Connect to Ensembl database
mart <- useEnsembl("ensembl", 
                   dataset = "hsapiens_gene_ensembl",
                   mirror = "useast")

#Fetch Hugo symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

#remove genes with no HUGO symbol
gene_map <- gene_map %>% filter(hgnc_symbol != "")

#Map HUGO names onto matrix
filtered_hugo <- filtered %>%
  rownames_to_column("ensembl_original") %>%
  mutate(ensembl_clean = str_remove(ensembl_original, "\\..*$")) %>%
  inner_join(gene_map, by = c("ensembl_clean" = "ensembl_gene_id")) %>%
  dplyr::select(-ensembl_original, -ensembl_clean) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  column_to_rownames("hgnc_symbol")

cat("Genes after HUGO mapping:", nrow(filtered_hugo), "\n")

#Log Transform
expr_matrix <- log2(filtered_hugo + 1)

#Filter to top 2000 most variable genes
gene_variance <- apply(expr_matrix, 1, var)
top_genes <- names(sort(gene_variance,, decreasing = TRUE))[1:2000]
expr_small <- expr_matrix[top_genes, ]



#saving as csv
write.csv(expr_small, "expression_matrix_hugo_filtered.csv")