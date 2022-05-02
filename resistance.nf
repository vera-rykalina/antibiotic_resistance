nextflow.enable.dsl = 2

process fastp {
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
  storeDir "${params.outdir}/fastp"
  input:
    path fastq

  output:
    path "${fastq.getSimpleName()}_fastp.fastq", emit: fastq_cut
    path "${fastq.getSimpleName()}_fastp.html", emit: fastpreport
  script:

    """
    fastp -i ${fastq} -o ${fastq.getSimpleName()}_fastp.fastq -h ${fastq.getSimpleName()}_fastp.html
    """
  }


  process fastqc {
    storeDir "${params.outdir}/fastqc"
    container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
    input:
      path fastq

    output:
      path "${fastq.getSimpleName()}_fastqc.html"
    script:
      """
      fastqc ${fastq}
      """
  }

process srst2 {
  storeDir "${params.outdir}/srst2"
  container "https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2"
  input:
    path combined
  output:
    path "${combined[0]}*.txt"
  script:
  """
  srst2 --input_se ${combined[0]} --output ${combined[0]} --gene_db ${combined[1]}
  """

}

process report {

  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path srst2
  output:
    path "final_report.txt"

  script:
  """
  cat *fullgenes__CARD_v3.0.8_SRST2__results.txt > final_report.txt
  """
}

workflow {
  fastqchannel = channel.fromPath("/home/cq/NGS/resistance_data/rawdata/*")
  fastpout = fastp(fastqchannel)
  fastqc_out = fastqc(fastpout.fastq_cut)
  reference = channel.fromPath("/home/cq/NGS/resistance_data/reference/*.fasta")
  reference.view()
  combined_channel = fastpout.fastq_cut.combine(reference)
  combined_channel.view()
  srst2_out = srst2(combined_channel)
  srst2_out.collect().view()
  final_rep = report(srst2_out.collect())

}
