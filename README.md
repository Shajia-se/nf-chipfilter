# nf-chipfilter

`nf-chipfilter` is a Nextflow DSL2 module for post-alignment BAM filtering in ChIP-seq workflows.

## What This Module Does

Filtering order:
1. Remove low-confidence / multi-mapping reads using MAPQ threshold (`samtools view -q`).
2. Optionally remove reads overlapping blacklist regions (`bedtools intersect -v`).
3. Remove mitochondrial reads (`chrM` / `MT`).

## Input

- Directory: `params.chipfilter_raw_bam`
- Optional: `params.samples_master` (CSV with `sample_id`, optional `enabled`) to restrict which samples are processed
- Input preference:
  - if `prefer_dedup=true` (default): use `*.dedup.bam`, fallback to `*.markdup.bam`
  - if `prefer_dedup=false`: use `*.markdup.bam`, fallback to `*.dedup.bam`

## Output

Under `${project_folder}/${chipfilter_output}`:
- `${sample}.nomulti.bam` + `.bai`
- `${sample}.noblack.bam` + `.bai` (if blacklist enabled)
- `${sample}.clean.bam` + `.bai`

## Key Parameters

- `chipfilter_raw_bam`: input BAM folder (usually `nf-picard/picard_output`)
- `samples_master`: optional sample whitelist source
- `chipfilter_output`: output folder name
- `prefer_dedup`: prefer dedup BAM as input (default: `true`)
- `mapq_threshold`: MAPQ filter cutoff (default: `4`)
- `blacklist_bed`: BED file for genomic blacklist filtering (optional)

## Run

```bash
nextflow run main.nf -profile hpc
```

```bash
nextflow run main.nf -profile hpc --mapq_threshold 4 --blacklist_bed /path/to/blacklist.bed
```

With sample restriction:

```bash
nextflow run main.nf -profile hpc \
  --chipfilter_raw_bam /path/to/nf-picard/picard_output \
  --samples_master /path/to/samples_master.csv \
  --mapq_threshold 4
```
