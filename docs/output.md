# czbiohub/sc2-ngs-analysis: Output

__DOCUMENTATION UNDER CONSTRUCTION__

This document describes the output produced by the pipeline. 

<!-- MarkdownTOC -->

- [Pipeline overview](#pipeline-overview)
- [Alignment to reference](#alignment-to-reference)
- [MultiQC](#multiqc)

<!-- /MarkdownTOC -->

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* Minimap2 - Align to reference
* Bcftools - Call variants
* [`assignclades.py`](../bin/assignclades.py) - Assign clades
* Bcftools - Search primer regions for variants
* BLAST - Find most similar sequences
* [`alignment_assembly_stats.py`](../bin/alignment_assembly_stats.py) - Per sample statistics
* Merge assemblies
* Merge stats
* Filter assemblies
* Rename assemblies
* Combine sequences
* Combine metadata
* Get CA sequences
* Make include file

## Alignment to reference

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
