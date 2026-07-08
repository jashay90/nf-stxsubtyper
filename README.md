# Nextflow Pipeline for Stx-subtyper
A nextflow pipeline to run multiple stx analyses on a genome and/or reads, and to generate a PDF summary report. This is still very much under development, so use at your own risk.

## Prerequisites
In order to run this pipeline, you must have [conda](https://docs.conda.io/projects/conda/en/stable/index.html) and [nextflow](https://www.nextflow.io/docs/latest/index.html) installed.

### Databases
You will also need to download the [CFIA/HC StxDB](https://github.com/OLC-Bioinformatics/STxDB), or another database which contains nucleotide sequence files for stx gene subunits, and nucleotide and protein sequence files for the full gene/protein. If this repository doesn't have prebuilt blast databases for the stx gene subunits, or a prebuilt kma database for the full genes, then you will have to make those.

To make the blast databases, use makeblastdb with default options.

To make the kma database, use kma_index with Kmersize 19.

## Running the pipeline
Go into the main directory for this repository, and run this command:

	nextflow run main.nf -c main.config

You may also run the pipeline from a different directory, in which case you would include the path to main.nf and main.config, but note that nextflow will create virtual environments for every directory in which you run this pipeline.

A number of parameters may be set either by changing the main.config file or by adding flags with two dashes to the command above.
### Parameters most likely to change
	genomes : assembly fasta files as input
	reads : paired read fastq files as input. If using, this should be formatted with quotes and a specific glob pattern for matching pairs, e.g. "path/to/*_{1,2}.fastq.gz"
	outdir: output directory
Note that in order to match reads with their associated assemblies, this pipeline is expecting reads filenames to have a prefix which matches the filename of the assembly, e.g. "name.fasta" matches with "name_{1,2}.fastq.gz"

### Database parameters
	blastdbprefix : the path and shared prefix to blast database files (e.g.. if your files are called nt_stx_1a_alleles.*, nt_stx_1b_alleles.*, nt_stx_2a_alleles.*, nt_stx_2b_alleles.*, then this parameter should be "stx_")
	blastdbsuffix : the shared suffix to blast database files (e.g.. if your files are called nt_stx_1a_alleles.*, nt_stx_1b_alleles.*, nt_stx_2a_alleles.*, nt_stx_2b_alleles.*, then this parameter should be "_alleles")
	refprefix : the shared prefix for the protein sequence files
	refsuffix : the shared suffix for the protein sequence files
	kmadb : the prefix for the kma database

