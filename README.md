# ChIP_seq

README for ChIP-seq pipeline

Designed and Created by Berkley Gryder and Hsien-Chao Chou at the National Cancer Institute

Purpose
ChIP-seq pipeline is a set of public and in-house scripts and programs connected together for processing ChIP-seq style data.

Structure
The pipeline consists of:
	1. DATA/ folders, where each experiment has its own folder for storing processed files (BAM files, TDF files, MACS peak calling output, HOMER motif output, etc)
	2. manage_sample/ folder, within which we maintain configuration files (chipseq.newBatch.conf) for global pipeline parameters and metadata for all samples in DATA
	3. scripts/ folder, within which we store all mature scripts for various parts of the pipeline
	4. projects/ folder, within which downstream cross-sample analysis (heatmaps, peak comparisons, plots, etc.) are stored in a folder for each project
	5. data_by_filetype/ folder, which contains reference data files used by the pipeline (bed files for blacklisted regions, for instance)

Requirements and Dependencies
The current pipeline is built to depend on biowulf2 on biowulf2.nih.gov;  all modules required (MACS2, ROSE, HOMER, etc.) are called from within various scripts underneath the runChipseqPipeline.pl
All that is additionally needed includes scripts (in scripts/ folder) and an indexed BAM file (DATA/Sample_LNCaP_H3K27ac_C_SRR2566837/Sample_LNCaP_H3K27ac_C_SRR2566837.bam)

Running the Pipeline
1. Assemble a sample file with metadata. example: ChIP_seq_samples_newBatch.txt

2. Assemble a configuration file (which references explicitly the metadata file to use). example: chipseq_newBatch.conf

3. Execute sbatch. example: sbatch -J newBatch  --partition=ccr --time=24:00:00 --export=config_file=/data/khanlab/projects/ChIP_seq/manage_samples/chipseq_newBatch.conf runChipseqPipeline.sh



