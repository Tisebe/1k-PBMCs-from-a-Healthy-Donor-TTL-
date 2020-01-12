# 1k-PBMCs-from-a-Healthy-Donor-TTL-
 This reproducible R Markdown analysis was created with workflowr (version 1.4.0). The Checks tab describes the reproducibility checks that were applied when the results were created. The Past versions tab lists the development history.  This vignette demonstrates how to manipulate bus format in R with BUSpaRse. The most recent version of bustools can generate gene count matrices from bus files more efficiently; the purpose of the separate implementation in BUSpaRse is for advanced users to experiment with new ways to collapse UMIs mapped to multiple genes and to adapt bus format to purposes other than single cell RNA-seq. This implementation is intended to facilitate exploration using R/Rcpp, which is easier to work with than C++.
 
 Setup

If you would like to rerun this notebook, you can git clone this repository or directly download this notebook from GitHub.
Install packages

We will be using the R packages below. BUSpaRse is not yet on CRAN or Bioconductor. For Mac users, see the installation note for BUSpaRse. The data, which is already in bus format, can be downloaded from the package TENxBUSData. Both TENxBUSData and BUSpaRse have been submitted to Bioconductor for review; the data for TENxBUSData can only be downloaded from Bioconductor devel (version 3.10), which requires R 3.6.0.

# Install devtools if it's not already installed
if (!require(devtools)) {
  install.packages("devtools")
}
# Install from GitHub
devtools::install_github("BUStools/BUSpaRse")
devtools::install_github("BUStools/TENxBUSData")

The package DropletUtils will be used to estimate the number of real cells as opposed to empty droplets. It’s on Bioconductor, and here is how it should be installed:

if (!require(BiocManager)) {
  install.packages("BiocManager")
}
# Install Bioconductor devel
BiocManager::install(version = "devel")
BiocManager::install("DropletUtils")

The other R packages below are on CRAN, and can be installed with install.packages.

library(BUSpaRse)
library(TENxBUSData)
library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(DropletUtils)
library(Matrix)
theme_set(theme_bw())

We will not get into the details of how to make the bus file with kallisto bus and bustools, as the data will be downloaded with TENxBUSData. TENxBUSData provides 5 different datasets, and we will use the PBMC 1k dataset here. The data from TENxBUSData contains the sorted bus file in text format. While the BUSpaRse package converts that text format bus file into gene count matrix, this text file can be loaded into R as a data frame for further exploration.

fn <- TENxBUSData("./output", dataset = "pbmc1k")

#> The dataset has already been downloaded. It is located in /Users/lambda/BUS_notebooks_R/output/out_pbmc1k

list.files(fn)

#> [1] "matrix.ec"         "output.sorted"     "output.sorted.txt"
#> [4] "transcripts.txt"

Explaining the output:

    matrix.ec: A text file with two columns. The first column is the 0 based index of equivalence classes. The second column is the set of transcripts (denoted by 0 based index based on order of appearance in the transcriptome fasta file) present in the corresponding equivalence class.
    output.sorted: The data represented in bus format, sorted by barcode, UMI, and equivalence class. This is a binary, so can’t be read into R with functions like read.table.
    output.sorted.txt: output.sorted converted into text format, so can be easily read into R for exploration.
    transcript.txt: A text file with one column, which is the transcripts present in the data, in the same order as in the transcriptome fasta file.

Sparse matrix
Map transcripts to genes

For the sparse matrix, most people are interested in how many UMIs per gene per cell, we here we will quantify this from the bus output, and to do so, we need to find which gene corresponds to each transcript. Remember in the output of kallisto bus, there’s the file transcripts.txt. Those are the transcripts in the transcriptome index. Information on which transcript corresponds to which gene can be directly retrieved from Ensembl.

tr2g <- transcript2gene(species = "Homo sapiens", 
                        kallisto_out_path = "./output/out_pbmc1k",
                        ensembl_version = 94)

#> Querying biomart for transcript and gene IDs of Homo sapiens

#> Sorting transcripts

head(tr2g)

#>           transcript              gene gene_name
#> 1: ENST00000632684.1 ENSG00000282431.1     TRBD1
#> 2: ENST00000434970.2 ENSG00000237235.2     TRDD2
#> 3: ENST00000448914.1 ENSG00000228985.1     TRDD3
#> 4: ENST00000415118.1 ENSG00000223997.1     TRDD1
#> 5: ENST00000631435.1 ENSG00000282253.1     TRBD1
#> 6: ENST00000390583.1 ENSG00000211923.1  IGHD3-10

Alternative ways of getting tr2g have been implemented in the BUSpaRse package. You may use tr2g_ensembl to query Ensembl with biomart to get transcript and gene IDs. If you use this method, then please make sure that the Ensembl version used in the query matches that of the transcriptome. This method is convenient for the user since you only need to input species names, but it can be slow since biomart database query can be slow. You may also use tr2g_gtf for GTF files and tr2g_gff3 for GFF3 files, which are more useful for non-model organisms absent from Ensemble. After calling the tr2g_* family of functions, you should sort the transcripts from those functions with sort_tr2g so the transcripts are in the same order as those in the kallisto index. Then the function save_tr2g_bustools can be used to save the tr2g data frame to a text file in the format required by bustools.
Make the sparse matrix

For 10x, we do have a file with all valid cell barcodes that comes with CellRanger. You need to install CellRanger to get this file, though you do not need to run CellRanger for this notebook. The whitelist is optional, so if you don’t have one, you may skip the whitelist step and the whitelist argument in the makr_sparse_matrix function.

# Copy v3 chemistry whitelist to working directory
cp ~/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt.gz \
./data/whitelist_v3.txt.gz

# Read in the whitelist
whitelist_v3 <- fread("./data/whitelist_v3.txt.gz", header = FALSE)$V1
length(whitelist_v3)

#> [1] 6794880

That’s an order of magnitude more than the 737K in v2 chemistry.

Now we have everything we need to make the sparse matrix. This function reads in output.sorted.txt line by line and processes them. It does not do barcode correction for now, so the barcode must exactly match those in the whitelist if one is provided. It took 5 to 6 minutes to construct the sparse matrix in the hgmm6k dataset, which has over 280 million lines in output.sorted.txt, which is over 9GB. Here the data set is smaller, and it takes less than a minute.

Note that the arguments est_ncells (estimated number of cells) and est_ngenes (estimated number of genes) are important. With the estimate, this function reserves memory for the data to be added into, reducing the need of reallocation, which will slow the function down. Since the vast majority of “cells” you get in this sparse matrix are empty droplets rather than cells, please put at least 200 times more “cells” than you actually expect in est_ncells.

If you do not have a whitelist of barcodes, then it’s fine; the whitelist argument is optional.

The function make_sparse_matrix can make the gene count matrix and the transcript compatibility count (TCC) matrix at the same time. For the purpose of this notebook, we only generate the gene count matrix. An upcoming notebook will demonstrate some more detailed analysis with a TCC matrix. See Ntranos et al. 2016 for more information about TCC matrices.

res_mat <- make_sparse_matrix("./output/out_pbmc1k/output.sorted.txt",
                              tr2g = tr2g, est_ncells = 3e5,
                              est_ngenes = nrow(tr2g),
                              whitelist = whitelist_v3, TCC = FALSE)

#> Reading matrix.ec
#> Processing ECs
#> Matching genes to ECs
#> Reading data
#> Read 5 million reads
#> Read 10 million reads
#> Read 15 million reads
#> Constructing gene count matrix

The matrix we get has genes in rows and barcode in columns. The row names are the gene IDs (not using human readable gene names since they’re not guaranteed to be unique), and the column names are cell barcodes.
Explore the data
Remove empty droplets

Cool, so now we have the sparse matrix. What does it look like?

dim(res_mat)

#> [1]  19821 216752

That’s way more cells than we expect, which is about 1000. So what’s going on?

How many UMIs per barcode?

tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>     1.00     1.00     2.00    46.49     6.00 54175.00

The vast majority of “cells” have only a few UMI detected. Those are empty droplets. 10x claims to have cell capture rate of up to 65%, but in practice, depending on how many cells are in fact loaded, the rate can be much lower. A commonly used method to estimate the number of empty droplets is barcode ranking knee and inflection points, as those are often assumed to represent transition between two components of a distribution. While more sophisticated method exist (e.g. see emptyDrops in DropletUtils), for simplicity, we will use the barcode ranking method here. However, whichever way we go, we don’t have the ground truth.

# Compute barcode rank
bc_rank <- barcodeRanks(res_mat)

qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")

The inflection point looks like a reasonable number of cells.

# Filter the matrix
res_mat <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
dim(res_mat)

#> [1] 19821  1169

Dimension reduction

seu <- CreateSeuratObject(res_mat, min.cells = 3) %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE)

See how number of total counts and number of genes expressed are distributed.

VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1)

Another QC plot

ggplot(seu@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "Total UMI counts per cell", y = "Number of genes detected")

seu <- RunPCA(seu, verbose = FALSE, npcs = 30)
ElbowPlot(seu, ndims = 30)

We can do Leiden clustering. Leiden is an improvement over Louvain that guarantees that clusters are well-connected on the k nearest neighbor graph.

# Leiden clustering
seu <- FindNeighbors(seu)

#> Computing nearest neighbor graph

#> Computing SNN

seu <- FindClusters(seu, algorithm = 4)

DimPlot(seu, reduction = "pca", pt.size = 0.5)

seu <- RunTSNE(seu, dims = 1:20, check_duplicates = FALSE)
DimPlot(seu, reduction = "tsne", pt.size = 0.5)



