# nf-chipfilter

A simple, portable BAM-filtering pipeline using Nextflow.
Designed to clean **markdup BAMs** (e.g., from nf-picard) for downstream ChIP-seq / ATAC-seq / RNA-seq analysis.

---

## 🎯 What this pipeline does

Runs **three standardized BAM-cleaning steps**:

### 1) Remove multimappers

Uses **MAPQ filtering** (default MAPQ ≥ 10):

```
samtools view -b -q 10
```

### 2) Remove blacklist regions

Uses **bedtools intersect -v** to exclude ENCODE blacklist peaks:

```
bedtools intersect -v -abam input.bam -b blacklist.bed
```

### 3) Remove mitochondrial reads

Keeps only non-mitochondrial chromosomes (e.g., chrM / MT removed):

```
samtools idxstats | grep -v "chrM"
samtools view -b input.bam <non-mito chr list>
```

Produces a final cleaned BAM ready for MACS2/3 peak calling.

---

## 🚀 Run on HPC (Slurm + Singularity)

```bash
nextflow run main.nf -profile hpc 
```

* Singularity image, e.g.:

```
singularity pull samtools-bedtools.sif docker://shajiase/samtools-bedtools:1.0
```

---

Each BAM passes through:

1. SAMPLE.nomulti.bam
2. SAMPLE.noblack.bam (only if a blacklist is provided)
3. SAMPLE.clean.bam ← **final cleaned BAM**

Final deliverables:

* `SAMPLE.clean.bam`
* `SAMPLE.clean.bam.bai`

Used for downstream peak calling (MACS2 / MACS3) or coverage tracks.

---

## 📂 Project structure

```
nf-chipfilter/
├── main.nf
├── nextflow.config
└── configs/
    ├── local.config
    └── slurm.config
```

---

## ✔️ Summary

This pipeline provides:

* Multi-mapper removal
* Blacklist removal
* Mitochondrial read removal
* Fully standardized clean BAM output
* Compatible with nf-bwa → nf-picard → nf-chipfilter → MACS2/3 analysis chain
