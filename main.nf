#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def chipfilter_output = params.chipfilter_output ?: "chipfilter_output"
def MAPQ             = params.mapq_threshold   ?: 4

process filter_multimappers {
  tag "${bam.simpleName}"
  stageInMode  'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${chipfilter_output}", mode: 'copy'

  input:
    path bam

  output:
    tuple path("${bam.simpleName}.nomulti.bam"),
          path("${bam.simpleName}.nomulti.bam.bai")

  script:
  """
  set -eux

  samtools view -b -q ${MAPQ} -o ${bam.simpleName}.nomulti.bam ${bam}
  samtools index ${bam.simpleName}.nomulti.bam
  """
}


process filter_blacklist {
  tag "${bam.simpleName}"
  stageInMode  'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${chipfilter_output}", mode: 'copy'

  input:
    tuple path(bam), path(bai) 

  output:
    tuple path("${bam.simpleName}.noblack.bam"),
          path("${bam.simpleName}.noblack.bam.bai")

  when:
    params.blacklist_bed

  script:
  """
  set -eux

  bedtools intersect -v -abam ${bam} -b ${params.blacklist_bed} > ${bam.simpleName}.noblack.bam
  samtools index ${bam.simpleName}.noblack.bam
  """
}

process remove_mito {
  tag "${bam.simpleName}"
  stageInMode  'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${chipfilter_output}", mode: 'copy'

  input:
    tuple path(bam), path(bai)

  output:
    tuple path("${bam.simpleName}.clean.bam"),
          path("${bam.simpleName}.clean.bam.bai")

  script:
  """
  set -eux

  keep_chrs=\$(samtools idxstats ${bam} | cut -f1 | egrep -v '^(chrM|MT|\\*)\$' | tr '\\n' ' ')

  samtools view -b ${bam} \$keep_chrs > ${bam.simpleName}.clean.bam
  samtools index ${bam.simpleName}.clean.bam
  """
}


workflow {

  def outdir = "${params.project_folder}/${chipfilter_output}"
  def prefer_dedup = (params.prefer_dedup == null) ? true : params.prefer_dedup
  def selectedSamples = null as Set
  def sampleMatches = { name, sid ->
    name == sid || name.startsWith("${sid}_")
  }

  def dedup_ch = Channel.fromPath("${params.chipfilter_raw_bam}/*.dedup.bam")
  def markdup_ch = Channel.fromPath("${params.chipfilter_raw_bam}/*.markdup.bam")
  def source_bam_ch = prefer_dedup ? dedup_ch.ifEmpty { markdup_ch } : markdup_ch.ifEmpty { dedup_ch }

  if (params.samples_master) {
    def master = file(params.samples_master)
    assert master.exists() : "samples_master not found: ${params.samples_master}"

    selectedSamples = [] as Set
    master.eachLine { line, n ->
      if (n == 1) return
      if (!line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (cols.size() < 1) return

      def sid = cols[0]
      if (!sid) return

      def enabled = cols.size() > 10 ? cols[10]?.toLowerCase() : ''
      if (enabled == '' || enabled == 'true') {
        selectedSamples << sid
      }
    }
    assert !selectedSamples.isEmpty() : "No enabled sample_id found in samples_master: ${params.samples_master}"
  }

  source_bam_ch
    .ifEmpty { exit 1, "ERROR: No input BAM found (*.dedup.bam or *.markdup.bam) under: ${params.chipfilter_raw_bam}" }
    .filter { bam ->
      if (selectedSamples == null) return true
      def sid = bam.simpleName.replaceFirst(/(?:\.sorted)?\.(dedup|markdup)$/, '')
      selectedSamples.any { wanted -> sampleMatches(sid, wanted as String) }
    }
    .filter { bam ->
      ! file("${outdir}/${bam.simpleName}.clean.bam.bai").exists()
    }
    .set { bam_ch }

  def nomulti_ch = filter_multimappers(bam_ch)

  def noblack_ch = params.blacklist_bed ? filter_blacklist(nomulti_ch) : nomulti_ch

  remove_mito(noblack_ch)
}
