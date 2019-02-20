# Classification comparison

Compare classifications of different tools in a metagenomics pipeline.  
Initial version: Compare [Contig Annotation Tool (CAT)](https://github.com/dutilh/CAT/)
to our in-house metagenomics pipeline, formerly known as "Pipeline Zonder Naam" (PZN).

-----

## Requirements

This project depends on:
 - [Snakemake](https://snakemake.readthedocs.io/en/stable/)  (tested with versions 5.2.2 and 5.3.0)
 - [conda](https://conda.io/en/latest/) (tested with version 4.6.4)
 - (optional) [DRMAA](https://pypi.org/project/drmaa/) (for use on a high-performance computing (HPC) cluster or grid computer)

Please make sure you have both installed before trying to run this software. 
Other software dependencies are installed with conda 
(for details see YAML files in `envs/`).

The pipeline has only been tested on **GNU/Linux**, specifically 
_Red Hat Enterprise Linux Server release 7.6 (Maipo)_. 
It is likely to work just as well on any conda and Snakemake compatible system.

The current configuration is made for use on a HPC cluster using 
[LSF](https://www.ibm.com/us-en/marketplace/hpc-workload-management) as 
scheduler with (Python) [DRMAA](https://pypi.org/project/drmaa/) as interface. 
If necessary, the configuration files may be edited to work on a different 
system.

-----

## User manual

_N.B. The manual will use conda for installing software. Other options exist,_
_but have not been tested._
_Also, this tool assumes you have already run the Jovian metagenomics pipeline_
_and have its results ready._

### 1. Initial setup

Install conda, following the instructions from 
https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation

Install Snakemake, for instance using a conda environment:  
`conda create -n snakemake snakemake=5.3.0 -c bioconda -c conda-forge`  
(version number may be omitted to install the latest available version)

Get the project code from GitHub:  
`git clone [link to project]`

### 2. Set the parameters to match your setup

_Optional: if not working on the same HPC cluster:_  
Edit the lines starting with "drmaa" and "jobname" in `conf/config.yaml` using 
a text editor. (Adjust to your own system, or remove the lines.  
Also, if you have more or fewer than 12 CPU cores on your machine, you can edit
the line "run_diamond_blastp: 12" under "threads:" accordingly.

_Optional: if you want to experiment with Prodigal/Diamond parameters:_
Edit the "Custom parameters" part (bottom) of `conf/parameters.yaml`.

Edit `conf/parameters.yaml` to point to the necessary files:
 - `source_dir` points to the Jovian folder
 - `work_dir` points to the folder where CAT's results and the comparison
 results will be stored 
 (for example: the folder where you downloaded this project)
 - `db_dir` points to the directory where the reference databases for CAT are
 stored 
 (for example: `/mnt/db/CAT`)
 - `samples` is a list of sample names (i.e. file prefixes) of the samples you
 want to analyse with CAT and compare to Jovian's Megablast (e.g. if you started
 with fastq files like `SRR7892426_R1.fastq` and `SRR7892426_R1.fastq`, then
 the sample name is "SRR7892426")
 
### 3. Run CAT and compare the results

(If using Snakemake from a conda environment, activate it.)  
`source activate snakemake`

Run snakemake:  
`snakemake --profile conf`

Wait until the analysis is finished and then browse the results in the
`results/` folder. Figures can be found in `results/figures/` and can be viewed
directly in a webbrowser. Tables in are in `results/tables/` and can for example
be used in small tests. E.g. the question "which virus species were identified
by both methods?" can be tested with 
`cut -f 1 results/tables/OVERALL.PZN-CAT.comparison.species.tsv | grep "virus"`.
(Such questions can be taylored to specific samples by changing "OVERALL" into
the sample name, to taxonomic ranks by changing "species" to the desired rank,
and to taxonomic names by changing "virus" to the desired name.)