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