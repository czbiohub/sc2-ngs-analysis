# czbiohub/sc2-ngs-analysis: Usage

## Table of contents
<!-- MarkdownTOC -->

- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Main arguments](#main-arguments)
  - [`-profile`](#-profile)
  - [`--sample_sequences`](#--sample_sequences)
  - [`--sample_metadata`](#--sample_metadata)
  - [`--blast_sequences`](#--blast_sequences)
  - [`--nextstrain_sequences`](#--nextstrain_sequences)
  - [`--nextstrain_metadata`](#--nextstrain_metadata)
  - [`--minLength`](#--minlength)
  - [`--maxNs`](#--maxns)
  - [`--clades`](#--clades)
  - [`--gisaid_names`](#--gisaid_names)
  - [`--qpcr_primers`](#--qpcr_primers)
  - [`--sample_vcfs`](#--sample_vcfs)
- [Reference genomes](#reference-genomes)
  - [`--ref`](#--ref)
  - [`--ref_gb`](#--ref_gb)
- [Job resources](#job-resources)
  - [Automatic resubmission](#automatic-resubmission)
  - [Custom resource requests](#custom-resource-requests)
- [AWS Batch specific parameters](#aws-batch-specific-parameters)
  - [`--awsqueue`](#--awsqueue)
  - [`--awsregion`](#--awsregion)
- [Other command line parameters](#other-command-line-parameters)
  - [`--outdir`](#--outdir)
  - [`--email`](#--email)
  - [`--email_on_fail`](#--email_on_fail)
  - [`-name`](#-name)
  - [`-resume`](#-resume)
  - [`-c`](#-c)
  - [`--custom_config_version`](#--custom_config_version)
  - [`--custom_config_base`](#--custom_config_base)
  - [`--max_memory`](#--max_memory)
  - [`--max_time`](#--max_time)
  - [`--max_cpus`](#--max_cpus)
  - [`--plaintext_email`](#--plaintext_email)
  - [`--monochrome_logs`](#--monochrome_logs)
  - [`--multiqc_config`](#--multiqc_config)

<!-- /MarkdownTOC -->


## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO nf-core: Document required command line parameters to run the pipeline-->

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run czbiohub/sc2-ngs-analysis -profile docker --sample_sequences '*.fa' --sample_metadeta sample_metadata.tsv --nextstrain_sequences sequences_yyyy-mm-dd.fasta --blast_sequences sequences_yy-mm-dd.fasta --nextstrain_metadata metadata_yyyy-mm-dd.tsv
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull czbiohub/sc2-ngs-analysis
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [czbiohub/sc2-ngs-analysis releases page](https://github.com/czbiohub/sc2-ngs-analysis/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
  * A configuration profile to run the pipeline using AWS Batch, specific to CZB AWS.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`czbiohub/sc2-msspe`](http://hub.docker.com/r/czbiohub/sc2-msspe/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`czbiohub/sc2-msspe`](http://hub.docker.com/r/czbiohub/sc2-msspe/)
* `test`
  * A profile with a complete configuration for automated testing
  * Does not include `--nextstrain_metadata` due to GISAID restrictions, so this file must be provided.

<!-- TODO nf-core: Document required command line parameters -->

### `--sample_sequences`
Use this to specify the location of your input SARS-CoV-2 sequence files. For example:

```bash
--sample_sequences 'path/to/data/sample_*.fa'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character

### `--sample_metadata`
Specify a file that contains the sample metadata, formatted to be compatible with Nextstrain. This file should specify sample collection dates.

```bash
--sample_metadata 'path/to/data/sample_metadata.tsv'
```

### `--blast_sequences`
Specify a file that contains sequences to build a BLAST database. This is normally the same file as `--nextstrain_sequences`.

```bash
--blast_sequences 'path/to/data/sequences.fasta'
```

### `--nextstrain_sequences`
Specify a file that contains the sequences that will be output into a file to be used as input for the `nextstrain/ncov` pipeline. This is normally the same file as `--blast_sequences`.

```bash
--nextstrain_sequences 'path/to/data/sequences.fasta'
```

### `--nextstrain_metadata`
Specify the nextstrain metadata file, normally retrieved from GISAID. This is used to filter for California sequences and to build a new metadata file that contains the sample metadata.

```bash
--nextstrain_metadata 'path/to/data/metadata.tsv'
```

### `--minLength`
Specify the minimum number of non-ambiguous/-missing/-gap characters needed for an assembly to pass filters.

```bash
--minLength 25000
```

By default, set to `25000` to conform to `nextstrain/ncov`. Note that an assembly must pass __both__ `--maxNs` and `--minLength` to pass filters.

### `--maxNs`
Specify the maximum number of N characters allowed for an assembly to pass filters.

```bash
--maxNs 6000
```

By default, set to `6000`. Note that an assembly must pass __both__ `--maxNs` and `--minLength` to pass filters.

### `--clades`
Specify the tab-separated clades file that contains clade definitions. This file needs to contains the columns `clade`, `gene`, `site`, `alt`. This is normally the same file as found in `nextstrain/ncov`.

```bash
--clades data/clades.tsv
```

By default, this is [`data/clades.tsv`](../data/clades.tsv)

### `--gisaid_names`
This should be a tab-separated file that contains the columns `sample_name`, `gisaid_name`. The file should map sample names, which correspond to the basename of the FASTA file, to GISAID names.

```bash
--gisaid_names gisaid_names.tsv
```

### `--qpcr_primers`
Specify a BED file with qPCR primer locations. This is needed to search for variants in these regions.

```bash
--qpcr_primers data/qpcr_primers.bed 
```

By default, this is [`data/qpcr_primers.bed`](../data/qpcr_primers.bed)

### `--sample_vcfs`
If the samples were run through `czbiohub/sc2-msspe-bioinfo/call_consensus.nf`, then each sample will have a corresponding existing VCF. This can be provided to prevent recomputation.

```bash
--sample_vcfs path/to/data/*.vcf
```

## Reference genomes


<!-- TODO nf-core: Describe reference path flags -->
### `--ref`
You can specify the full path to your reference genome when you run the pipeline:

```bash
--ref '[path to Fasta reference]'
```

By default, the pipeline uses Wuhan-Hu-1, found in the file [`data/MN908947.3.fa`](../data/MN908947.3.fa).

###  `--ref_gb`
A Genbank reference file is needed to assign clades.

```bash
--ref_gb '[path to Genbank reference]'
```

By default, the pipeline uses Wuhan-Hu-1, found in the file [`data/MN908947.3.gb`](../data/MN908947.3.gb)

## Job resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `us-west-2` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`
This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`
Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`
If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`
Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.
