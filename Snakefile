from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

configfile: "config.yaml"
# SAMPLES=["SRR1543136", "SRR1543137"]

# All rule
rule all:
	input:
		r_image="results/image.Rdata",
		dendogram="results/dendogram.pdf",
		pca="results/pca.pdf",
		dge_table="results/topgenes.csv",
		ma_plot="results/ma_plot.pdf"

rule fasterq_dump:
	output:
		"reads/raw/{sample}.fastq"
	shell:
		"fasterq-dump -o {output} {wildcards.sample}"

# Download reference Fasta from Ensembl
rule reference_fasta:
	input:
		FTP.remote(config["reference"]["fasta"])
	output:
		"reference/fasta/genome.fa"
	shell:
		"gzip -d -c {input} > {output}"

# Download reference GTF from Ensembl
rule reference_gtf:
	input:
		FTP.remote(config["reference"]["gtf"])
	output:
		"reference/gtf/genome.gtf"
	shell:
		"gzip -d -c {input} > {output}"

# Generate genome indexes
rule genome_index:
	input:
		fasta="reference/fasta/genome.fa",
		gtf="reference/gtf/genome.gtf"
	output:
		"reference/index/chrLength.txt",
		"reference/index/chrNameLength.txt",
		"reference/index/chrName.txt",
		"reference/index/chrStart.txt",
		"reference/index/exonGeTrInfo.tab",
		"reference/index/exonInfo.tab",
		"reference/index/geneInfo.tab",
		"reference/index/Genome",
		"reference/index/genomeParameters.txt",
		"reference/index/Log.out",
		"reference/index/SA",
		"reference/index/SAindex",
		"reference/index/sjdbInfo.txt",
		"reference/index/sjdbList.fromGTF.out.tab",
		"reference/index/sjdbList.out.tab",
		"reference/index/transcriptInfo.tab"
	threads: 4
	shell:
		"STAR --runThreadN {threads} "
		"--runMode genomeGenerate "
		"--genomeDir reference/index "
		"--genomeFastaFiles {input.fasta} "
		"--sjdbGTFfile {input.gtf} "
		"--sjdbOverhang 75 "
		"--genomeSAindexNbases 12"

# Perform alignment
rule mapped_reads:
	input:
		"reference/index/chrLength.txt",
		"reference/index/chrNameLength.txt",
		"reference/index/chrName.txt",
		"reference/index/chrStart.txt",
		"reference/index/exonGeTrInfo.tab",
		"reference/index/exonInfo.tab",
		"reference/index/geneInfo.tab",
		"reference/index/Genome",
		"reference/index/genomeParameters.txt",
		"reference/index/Log.out",
		"reference/index/SA",
		"reference/index/SAindex",
		"reference/index/sjdbInfo.txt",
		"reference/index/sjdbList.fromGTF.out.tab",
		"reference/index/sjdbList.out.tab",
		"reference/index/transcriptInfo.tab",
		reads="reads/raw/{sample}.fastq"
	threads: 4
	output:
		"reads/mapped/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	shell:
		"STAR --runThreadN {threads} "
		"--readFilesIn reads/raw/{wildcards.sample}.fastq "
		"--genomeDir reference/index "
		"--outFileNamePrefix reads/mapped/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate "
		"--outSAMunmapped Within "
		"--outSAMattributes Standard "
		"--alignIntronMax 20000"
		
rule exp_analysis:
	input:
		samples_table=config["samples-table"],
		bam_files=expand("reads/mapped/{sample}/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"]),
		gtf="reference/gtf/genome.gtf"
	output:
		r_image="results/image.Rdata",
		dendogram="results/dendogram.pdf",
		pca="results/pca.pdf",
		dge_table="results/topgenes.csv",
		ma_plot="results/ma_plot.pdf"
	script:
		"scripts/exp-analysis.R"	