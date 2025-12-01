#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def chipfilter_output = params.chipfilter_output ?: "chipfilter_output"

def MAPQ = params.mapq_threshold ?: 10


process filter_multimappers {
  tag "${bam.simpleName}"
  stageInMode  'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${chipfilter_output}", mode: 'copy'

  input:
    path bam

  output:
    path "${bam.simpleName}.nomulti.bam"
    path "${bam.simpleName}.nomulti.bam.bai"

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
    path bam

  output:
    path "${bam.simpleName}.noblack.bam"
    path "${bam.simpleName}.noblack.bam.bai"

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
    path bam

  output:
    path "${bam.simpleName}.clean.bam"
    path "${bam.simpleName}.clean.bam.bai"

  script:
  """
  set -eux

  keep_chrs=\$(samtools idxstats ${bam} | cut -f1 | egrep -v '^(chrM|MT|\*)\$' | tr '\\n' ' ')
  samtools view -b ${bam} \$keep_chrs > ${bam.simpleName}.clean.bam

  samtools index ${bam.simpleName}.clean.bam
  """
}


workflow {

  def outdir = "${params.project_folder}/${chipfilter_output}"

  def bams = Channel
    .fromPath("${params.chipfilter_raw_bam}/*.markdup.bam")
    .filter { bam ->
      ! file("${outdir}/${bam.simpleName}.clean.bam.bai").exists()
    }

  def nomulti = filter_multimappers(bams)

  def noblack = nomulti
  if( params.blacklist_bed ) {
    noblack = filter_blacklist(nomulti)
  }

  remove_mito(noblack)
}
