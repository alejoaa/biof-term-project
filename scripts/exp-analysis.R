#!/usr/bin/env R

# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

if (!requireNamespace("GenomicFeatures", quietly = TRUE))
BiocManager::install("GenomicFeatures")

if (!requireNamespace("Rsamtools", quietly = TRUE))
BiocManager::install("Rsamtools")

if (!requireNamespace("GenomicAlignments", quietly = TRUE))
BiocManager::install("GenomicAlignments")

if (!requireNamespace("BiocParallel", quietly = TRUE))
BiocManager::install("BiocParallel")

if (!requireNamespace("rtracklayer", quietly = TRUE))
BiocManager::install("rtracklayer")

if (!requireNamespace("DESeq2", quietly = TRUE))
BiocManager::install("DESeq2")

if (!requireNamespace("gplots", quietly = TRUE))
install.packages("gplots")

# Import packages
suppressMessages(library("Rsamtools"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("GenomicAlignments"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("DESeq2"))
suppressMessages(library("gplots"))
suppressMessages(library("rtracklayer"))

# Set up features for SummarizeOverlaps
gff0 <- import(snakemake@input[["gtf"]])

idx <- mcols(gff0)$source == "protein_coding" & mcols(gff0)$type == "exon" & (seqnames(gff0) %in% c("X",  "2L", "2R", "3L", "3R", "4", "YHet"))

gff <- gff0[idx]

genes <- split(gff, mcols(gff)$gene_id)

# Import table with columns FileName SRRNumber GEOAccession Status SampleName
cat("- Importing samples table\n")
samples_table <- read.table(file=snakemake@input[["samples_table"]], sep="\t", header=TRUE)

# Create txdb object
# cat("- Creating txdb object\n")
# txdb <- makeTxDbFromGFF(file=snakemake@input[["gtf"]], format="gtf", organism="Drosophila melanogaster")

# organize exons by gene
# cat("- Organizing exons by gene\n")
# eByg <- exonsBy(txdb, by="gene")

# Import bam files
cat("- Importing bam files\n")
cat(snakemake@input[["bam_files"]])
# bam_names <- as.character(samples_table$FileName)
bam_files <- BamFileList(snakemake@input[["bam_files"]])

# SummarizeOverlaps
cat("- Running SummarizeOverlaps\n")
register(MulticoreParam(multicoreWorkers()))

se <- summarizeOverlaps(features=genes, reads=bam_files, mode="Union", singleEnd=TRUE, ignore.strand=FALSE)

# Format se
colData(se) <- DataFrame(samples_table)
colnames(se) <- samples_table$Accession

# rowData(se) <- names(rowRanges(se))
# Alternative to the above method
# colnames(rowRanges(se)) <- "id" #the colname is changed to "FBtr"

#---------------------------------
# Differential expression analysis
#---------------------------------

cat("- Performing DGE\n")
# Build DESeq Data set
dds <- DESeqDataSet(se, design= ~ Enviroment)

# Regularized log ratio
rld <- rlog(dds)

# Measure euclidean distance
d <- dist(t(assay(rld)))

# Agglomerate using complete distance
hc <- hclust(d)

# Plot dendogram
dend = as.dendrogram(hc)
pdf(snakemake@output[["dendogram"]])
plot(dend)
dev.off()

# PCA
pdf(snakemake@output[["pca"]])
plotPCA(rld, intgroup = "Enviroment")
dev.off()

# Generate results table
dds <- DESeq(dds)
resD <- results(dds, alpha=0.05)

# Save DE genes in csv file
resDSort <- resD[order(resD$padj),]
topDESeq2 <- resDSort[1:395,]
write.csv(topDESeq2, file=snakemake@output[["dge_table"]])

# MA plot
pdf(file=snakemake@output[["ma_plot"]])
plotMA(resD, ylim=c(-7,7))
dev.off()

cat("- DGE done\n")
save.image(file=snakemake@output[["r_image"]])