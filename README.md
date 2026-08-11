# nf-chipfilter

`nf-chipfilter` filters Picard-processed BAM files for ChIP-seq peak calling.

This module is normally the fifth step in the ChIP-seq pipeline:

```text
nf-fastqc -> nf-fastp -> nf-bwa -> nf-picard -> nf-chipfilter -> nf-macs3
```

`nf-picard` produces deduplicated or duplicate-marked BAM files. `nf-chipfilter` removes low-confidence alignments by MAPQ and reports mitochondrial read burden.

## What It Does

1. Finds Picard output BAM files.
2. Prefers `*.dedup.bam` by default, with fallback to `*.markdup.bam`.
3. Optionally restricts samples using `samples_master`.
4. Filters reads using `samtools view -q`.
5. Indexes the filtered BAM.
6. Writes a small QC table with MAPQ-filtered reads and mitochondrial read estimates.

## Before You Run

You need:

- Nextflow installed
- BAM files from `nf-picard`
- For local runs: Docker available
- For HPC runs: Slurm and Singularity available
- For HPC runs: the Samtools/Bedtools Singularity image path in `configs/slurm.config` must exist

Default HPC notification email:

```text
molendo.hpc@gmail.com
```

## Input

The input folder should usually be the `nf-picard` output folder.

Expected input files:

```text
sample.sorted.dedup.bam
sample.sorted.dedup.bam.bai
```

or, if duplicate removal was disabled in Picard:

```text
sample.sorted.markdup.bam
sample.sorted.markdup.bam.bai
```

The module reads:

```bash
--chipfilter_raw_bam /path/to/picard_output
```

Optional sample restriction:

```bash
--samples_master /path/to/samples_master.csv
```

When `samples_master` is provided:

- It must contain a `sample_id` column.
- If it contains an `enabled` column, rows with `enabled=false` are skipped.
- Empty `enabled` values are treated as enabled.

## Input Preference

Default:

```bash
--prefer_dedup true
```

This uses `*.dedup.bam` when available and falls back to `*.markdup.bam`.

Alternative:

```bash
--prefer_dedup false
```

This uses `*.markdup.bam` when available and falls back to `*.dedup.bam`.

For this ChIP-seq pipeline, deduplicated BAM is usually preferred before peak calling.

## Output

Results are written to:

```text
${project_folder}/${chipfilter_output}/
```

Example:

```text
/path/to/output_project/chipfilter_output/
  WT_rep1.sorted.dedup.nomulti.bam
  WT_rep1.sorted.dedup.nomulti.bam.bai
  WT_rep1.sorted.dedup.chipfilter.stats.tsv
```

Output files:

- `*.nomulti.bam`: MAPQ-filtered BAM
- `*.nomulti.bam.bai`: BAM index
- `*.chipfilter.stats.tsv`: compact filtering/QC summary

The stats table contains:

```text
sample_id
mapq_threshold
nomulti_reads
mito_reads
pct_mito
clean_reads_estimated
```

## Recommended HPC Run

From inside the `nf-chipfilter` folder:

```bash
cd /path/to/nf-chipfilter

nextflow run main.nf -profile hpc \
  --chipfilter_raw_bam /path/to/picard_output \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --chipfilter_output chipfilter_output \
  --mapq_threshold 24
```

Resume a previous run:

```bash
nextflow run main.nf -profile hpc -resume \
  --chipfilter_raw_bam /path/to/picard_output \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --chipfilter_output chipfilter_output \
  --mapq_threshold 24
```

Override the HPC notification email:

```bash
nextflow run main.nf -profile hpc \
  --chipfilter_raw_bam /path/to/picard_output \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --mail_user molendo.hpc@gmail.com
```

## Local Test Run

Use local mode only for small test BAM files:

```bash
nextflow run main.nf -profile local \
  --chipfilter_raw_bam /path/to/test_picard_output \
  --project_folder ./test_output
```

## Key Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `chipfilter_raw_bam` | Input folder containing Picard BAM files | `null` |
| `samples_master` | Optional CSV used to select enabled sample IDs | `null` |
| `chipfilter_output` | chipfilter output subfolder | `chipfilter_output` |
| `prefer_dedup` | Prefer `*.dedup.bam` over `*.markdup.bam` | `true` |
| `mapq_threshold` | MAPQ cutoff for keeping reads | `24` |
| `project_folder` | Base output folder | current folder |
| `cpus` | CPUs per default task | `2` |
| `memory` | Memory per default task | `8 GB` |
| `time` | Runtime limit per default task | `8h` |
| `mail_user` | HPC notification email | `molendo.hpc@gmail.com` |

## Existing Results Are Skipped

For each input BAM, the module checks whether both expected outputs already exist:

```text
sample.nomulti.bam
sample.nomulti.bam.bai
sample.chipfilter.stats.tsv
```

If all expected files exist, that sample is skipped. If any file is missing, filtering runs again for that sample.

## How To Check Results

Start with:

```text
${project_folder}/${chipfilter_output}/sample.chipfilter.stats.tsv
```

Important values:

- `nomulti_reads`: reads retained after MAPQ filtering
- `mito_reads`: mitochondrial reads after MAPQ filtering
- `pct_mito`: mitochondrial percentage
- `clean_reads_estimated`: approximate usable non-mitochondrial reads

The `*.nomulti.bam` and `*.nomulti.bam.bai` files are the inputs for the next module, usually `nf-macs3`.

## Troubleshooting

If the run fails:

1. Check the main Nextflow log:

```bash
less .nextflow.log
```

2. Check the failed task error file:

```bash
less work/<hash>/.command.err
```

3. Common problems:

- `--chipfilter_raw_bam` does not point to the `nf-picard` output folder.
- No `*.dedup.bam` or `*.markdup.bam` files are present.
- `sample_id` in `samples_master` does not match the BAM prefix.
- Chromosome names do not use `chrM` or `MT`; mitochondrial counts may be zero.
- `configs/slurm.config` points to a missing Singularity image.
- The HPC bind path in `extra_mounts` does not include the BAM or output location.
- Docker is not running for local mode.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
