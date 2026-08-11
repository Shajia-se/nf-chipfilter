#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def chipfilter_output = params.chipfilter_output ?: "chipfilter_output"
def MAPQ             = params.mapq_threshold   ?: 24

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


process mito_stats {
  tag "${bam.simpleName}"
  stageInMode  'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${chipfilter_output}", mode: 'copy'

  input:
    tuple path(bam), path(bai)

  output:
    path("*.chipfilter.stats.tsv")

  script:
  def sample_base = bam.simpleName.replaceFirst(/\.nomulti$/, '')
  """
  set -eux

  nomulti_reads=\$(samtools view -c -F 260 ${bam})
  mito_reads=\$(samtools idxstats ${bam} | awk '\$1 == "chrM" || \$1 == "MT" {sum += \$3} END{print sum+0}')
  clean_reads_estimated=\$(( nomulti_reads - mito_reads ))

  pct_mito=NA
  if [[ "\$nomulti_reads" -gt 0 ]]; then
    pct_mito=\$(awk -v a="\$mito_reads" -v b="\$nomulti_reads" 'BEGIN{printf "%.2f", (a/b)*100}')
  fi

  cat > ${sample_base}.chipfilter.stats.tsv << TSV
sample_id	mapq_threshold	nomulti_reads	mito_reads	pct_mito	clean_reads_estimated
${sample_base}	${MAPQ}	\$nomulti_reads	\$mito_reads	\$pct_mito	\$clean_reads_estimated
TSV
  """
}


workflow {

  if (!params.chipfilter_raw_bam) {
    exit 1, "ERROR: --chipfilter_raw_bam must be provided. This should usually be the nf-picard output folder."
  }

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
    master
      .readLines()
      .findAll { it?.trim() }
      .with { lines ->
        def header = lines[0].split(',', -1)*.trim()
        def sampleIdx = header.indexOf('sample_id')
        def enabledIdx = header.indexOf('enabled')
        assert sampleIdx >= 0 : "samples_master must contain a sample_id column: ${params.samples_master}"

        lines.drop(1).each { line ->
          def cols = line.split(',', -1)*.trim()
          def sid = cols.size() > sampleIdx ? cols[sampleIdx] : ''
          if (sid) {
            def enabled = enabledIdx >= 0 && cols.size() > enabledIdx ? cols[enabledIdx]?.toLowerCase() : ''
            if (enabled == '' || enabled == 'true') {
              selectedSamples << sid
            }
          }
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
      !file("${outdir}/${bam.simpleName}.nomulti.bam").exists() ||
      !file("${outdir}/${bam.simpleName}.nomulti.bam.bai").exists() ||
      !file("${outdir}/${bam.simpleName.replaceFirst(/\\.nomulti$/, '')}.chipfilter.stats.tsv").exists()
    }
    .set { bam_ch }

  def nomulti_ch = filter_multimappers(bam_ch)
  mito_stats(nomulti_ch)
}
