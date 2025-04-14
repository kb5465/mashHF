#!/usr/bin/env Rscript

############################################################
# mashHF: Filter REGENIE Set Inclusion File by Gene Chunk
# Author: Kiran Biddinger (adapted from Sean Jurgens)
# Source: Adapted from 
#   https://github.com/seanjosephjurgens/TOPMed_AFib_pipeline/blob/patch-1/DNANexus/REGENIE_pipeline/v2/filter_setinclusion_files.R
#
# Description:
#   This script filters a REGENIE step 2 set inclusion file based on:
#     1. A user-specified list of target genes
#     2. Optional restriction to specific tissues (e.g., canonical, LV)
#     3. Optional filtering by chunk index for parallel processing
#
#   The output is a filtered set list used in gene-level association testing.
############################################################

# ----------------------------
# 1. Parse Input Arguments
# ----------------------------
args <- commandArgs(TRUE)
regenie_setinclusionfile      <- as.character(args[1])
gene_chunk_num                <- as.numeric(args[2])
gene_chunk_total              <- as.numeric(args[3])
genes_to_keep                 <- as.character(args[4])
tissues_to_keep               <- as.character(args[5])
regenie_setinclusionfile_out  <- as.character(args[6])

# ----------------------------
# 2. Load Required Libraries and Paths
# ----------------------------
library(data.table)
.libPaths(c("rpackages4_1_3", .libPaths()))  # Add custom library path if needed

# ----------------------------
# 3. Load Set Inclusion File
# ----------------------------
group <- fread(regenie_setinclusionfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)

# Extract ENSG IDs (before "__") from each grouping line
genes <- unique(gsub("__.*", "", group$V1))

# ----------------------------
# 4. Subset to Specified Genes
# ----------------------------
genes_inclusion <- strsplit(genes_to_keep, ",")[[1]]
genes <- genes[genes %in% genes_inclusion]

# ----------------------------
# 5. Chunk the Genes and Select One Chunk
# ----------------------------
chunk_size <- ceiling(length(genes) / gene_chunk_total)
splitz <- split(seq_along(genes), ceiling(seq_along(genes) / chunk_size))
genes_chunk <- genes[splitz[[gene_chunk_num]]]

# ----------------------------
# 6. Filter Group File by Chunk Genes
# ----------------------------
group <- group[gsub("__.*", "", group$V1) %in% genes_chunk, , drop = FALSE]
colnames(group) <- "V1"

# ----------------------------
# 7. Filter by Tissues of Interest
# ----------------------------
tissues_to_keep_vec <- strsplit(tissues_to_keep, ",")[[1]]
group <- group[gsub(".*__", "", group$V1) %in% tissues_to_keep_vec, , drop = FALSE]

# ----------------------------
# 8. Sort and Write Output
# ----------------------------
group <- group[order(group$V1), , drop = FALSE]
write.table(group, file = regenie_setinclusionfile_out, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
