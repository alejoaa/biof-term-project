# Differential Gene Expression Analysis Workflow

# Background and rationale
The widespread use of RNA-Seq in recent years has generated a large amount of open-source datasets of transcriptomes from a wide variety of organisms, including laboratory model species and human. One of the most common analyses with RNA-Seq data is the identification of differentially expressed genes between 2 or more groups, conditions, or treatments. For example, differential gene expression (DGE) analysis is particularly useful for assessing the effect of [chemotherapies in cancer](https://doi.org/10.1093/nar/gkw797), by studying changes in the transcriptome between the responsive and non-responsive treatment groups. Another application is when studying [model organisms](https://doi.org/10.1186/s12864-015-2353-z), by assessing the background variability of the transcriptome for individuals of different sex or raised in distinct environments. The results from a DGE analysis come in the form of a table of genes with adjusted p-values that indicate their statistical significance across the groups, and the logfold change which measures the amount and direction of the expression changes (_i.e._ if a given gene is up or downregulated). Other results can include visualizations, such as an MA plot or heatmap.

The general bioinformatic steps for performing DGE are:
1. Generate BAM alignments
2. Generate counts matrix for genes and samples
3. Perform differential gene expression analysis on the matrix

Each step can make use of different tools developed by the open-source community. For instance, alignments can be created with TopHat2 or STAR, counts obtained with SummarizeOverlaps or FeatureCounts, and DGE can be assessed with DESeq2 or edgeR. While this diversity in tools presents an opportunity for fine-tuning DGE analysis in terms of algorithms, statistical models, and computing resources, it also results in a barrier for reproducibility as the combination of tools, versions, and operating systems makes it difficult to run a pipeline as is in various computing contexts. As a consequence, the validation of results that comes from DGE analysis is hampered, as well as the building of a generic workflow that is robust across operating systems and package versions. This is particularly pertinent given that different strategies for analyzing RNA-Seq data have been shown to perform with different [rediscovery rates (RDRs) for top-ranked genes](https://doi.org/10.3389/fgene.2019.01331). Additionally, many DGE studies do not provide information about the specific commands used for the tools employed, nor the parameters for the options that each of them have.

In order to approach a solution for this issue, the current work defines a workflow for DGE analysis using Snakemake, the STAR aligner, SummarizeOverlaps, and DESeq2. Most packages are installed through `conda` or automatically in R, so the only technical requirement is `conda V4.8.4` or higher (earlier version may work as well although they were not tested). Installation instructions for `conda` are provided by their official documentation for [macOS, Linux, and Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation). It is important to note that the present workflow has been tested in CentOS 7 and macOS10.15.7 (Catalina). Running it on different OS versions may require small tweaks, but the commands should work for Unix systems in general without issues. The main dependencies, for which instructions will be provided in the usage section, are the following:

* `r-base =4.0.3`
* `pip =20.2.4`
* `snakemake = 5.30.1`
* `sra-tools =2.8.0`
* `star =2.7.6`

An overview of the workflow can be represented with a DAG. The following DAG contains only 4 of the 12 samples used in the example run for simplicity.

![alt text](https://github.com/alejoaa/biof-term-project/blob/main/dag.png)

Therefore, the main steps of the workflow are:
1. Download reference files: FASTA and GTF
2. Download reads from SRA with `fasterq-dump`
3. Build genome index with `STAR`
4. Map reads to reference using `STAR` and generate BAM files
5. Perform DGE analysis with custom `exp-analysis.R` script, which uses SummarizeOverlaps and DESeq2

# Usage
First, clone the git repository to your local machine and cd into the folder.
```
# Clone git repo
git clone https://github.com/alejoaa/biof-term-project.git
cd biof-term-project
```

Next, create an environment from the `environment.yaml` file using `conda`. In this example, the created env is called `exp-analysis`, but you can replace with a preferred name.
```
# Create conda environment
conda env create --name exp-analysis -f enviroment.yaml
conda activate exp-analysis
```

Given that some incompatibilities where found between the requirements of `R 4.0` and `Snakemake`, we’ll use `pip` for installing the last one.  `pip` and `R` were installed by `conda` in the previous step.
```
# Install Snakemake
pip install snakemake
# Install ftputil required by Snakemake later
pip install ftputil
```

With the previous steps all requirements have been met for running the workflow. Finally, we run `Snakemake` with the number of cores desired. For most personal computers these are between 2 - 4. Also, note that the task used for the run is called `all`. This follows the best practice of providing a single task that generates all of the outputs for the workflow.
```
# Run Snakemake with a specific number of cores
# Option -p prints the commands that are executed
snakemake --cores 4 -p all
```

# Input
The input data for the workflow consists of :
* 12 fastq files of /Drosophila melanogaster/ transcriptomes
* TSV table containing the metadata for the samples
* Reference genome fasta and GTF
* `config.yaml` file

The reads  used in the present run where generated by [Lin, Y., Golovnina, K., Chen, ZX. _et al._](https://doi.org/10.1186/s12864-015-2353-z) as part of a study of 726 _D. melanogaster_ individuals. Their project aimed to compare the differences in the transcriptomes across 3 factors: genotype, environment, and sex. The complete dataset can be found under the NCBI Gene Expression Omnibus (GEO) [GSE60312](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi). The 12 samples used for the workflow belong to the same genotype (_Drosophila_ Genetic Reference Panel DGRP-804) and sex (male), while differing in the environments that the organisms were raised in (E1 or E2):

|Accession|Genotype|Sex|Environment|
|---|---|---|---|
|SRR1543294|DGRP-804|M|E1|
|SRR1543295|DGRP-804|M|E1|
|SRR1543323|DGRP-804|M|E1|
|SRR1543324|DGRP-804|M|E1|
|SRR1543388|DGRP-804|M|E1|
|SRR1543389|DGRP-804|M|E1|
|SRR1543509|DGRP-804|M|E2|
|SRR1543514|DGRP-804|M|E2|
|SRR1543523|DGRP-804|M|E2|
|SRR1543537|DGRP-804|M|E2|
|SRR1543564|DGRP-804|M|E2|
|SRR1543580|DGRP-804|M|E2|

The TSV table corresponds to the `samples-table.tsv` file in the repository, and it has the same structure as the table above. The reference fasta and GTF files are specified as an FTP link to the Ensembl database in the `config.yaml` file. This file also indicates the list of SRA accessions, which must be in the same order as the rows in `samples-table.tsv`.  The reference files as well as fastq reads are automatically downloaded by the workflow into the folders `reference/fasta/`, `reference/gtf`, and `reads/raw`. The workflow is not restricted to these 12 samples, and if the user wants to include more from the study, they can do so by adding the samples and relevant metadata in `config.yaml` and `samples-table.tsv`.

# Output
The outputs of the workflow can be found under the `results` folder and are generated during the execution of the last task in the workflow (`exp_analysis`). This task runs the `exp-analysis.R` script under the `scripts` folder, which takes as input the mapped reads in `BAM` format and the reference `GTF` file for performing the DGE analysis. The expected results can be found in the `example-results` folder. In brief, the results consist of:
* Dendogram: `dendogram.pdf`
* PCA plot: `pca.pdf`
* Table with top differentially expressed genes: `topgenes.cvs`
* MA plot: `ma_plot.pdf`
* R image: `image.Rdata`

The *dendogram* and *PCA* plots are based on the regularized-logarithm transformation of the raw counts of reads per gene. The dendogram uses hierarchical clustering with Euclidean distance. For the samples used in the workflow, it is possible to note that the dendogram splits the ones in E1 and E2 into 2 different clusters. The same discrimination is observed in the PCA plot, where the main principal component (PC) which accounts for 47% of the variance can explained by the differences of conditions of E1 and E2. Interestingly, the individuals in E1 display more variability between them than the ones in E1, according to the second PC.

![alt text](https://github.com/alejoaa/biof-term-project/blob/main/example-results/dendogram.png)
![alt text](https://github.com/alejoaa/biof-term-project/blob/main/example-results/pca.png)

The table of differentially expressed genes contains the top 395 ordered by their p-adjusted value. The table also contains the following columns:
* `baseMean`: base mean over all genes
* `log2FoldChange`: log2 fold change of E1 vs E2 (E1 is the reference)
* `lfcSE`: standard error of the log2 fold change
* `stat`: Wald statistic of E1 vs E2
* `pvalue`: Wald test p-value of E1 vs E2
* `padj`: adjusted p-values

The top 5 differentially expressed genes by their FlyBase IDs are _FBgn0031273_, _FBgn0004551_, _FBgn0086783_, _FBgn0040519_, _FBgn0038181_.

The MA plot displays the genes as dots where the y-axis is the log fold change and the x-axis is the mean of normalized counts. The dots that have blue fill are the genes that are a statistically significant differential expression between E1 and E2.
![alt text](https://github.com/alejoaa/biof-term-project/blob/main/example-results/ma_plot.png)

Finally, the last result generated by the workflow is the R image, which contains all the variables used during the run of `exp-analysis.R`. The main variables are:
* `se`: output of SummarizeOverlaps
* `rld`: regularized-logarithm transformed counts
* `resD`: result of DESeq2 after performing DGE analysis
