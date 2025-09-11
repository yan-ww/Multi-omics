version 1.0

workflow scRNAseq_Pipeline {
  input {
    Array[File] fastq_files_1
    Array[File] fastq_files_2
    File reference_genome
    File gtf_annotation
    String sample_prefix
  }

  call FastQC {
    input:
      fastq_files_1 = fastq_files_1,
      fastq_files_2 = fastq_files_2
  }

  call STAR_Alignment {
    input:
      fastq_files_1 = fastq_files_1,
      fastq_files_2 = fastq_files_2,
      reference_genome = reference_genome,
      gtf_annotation = gtf_annotation,
      sample_prefix = sample_prefix
  }

  call FeatureCounts {
    input:
      bam_files = STAR_Alignment.bam_files,
      gtf_annotation = gtf_annotation,
      sample_prefix = sample_prefix
  }

  output {
    Array[File] fastqc_reports = FastQC.reports
    Array[File] aligned_bams = STAR_Alignment.bam_files
    File count_matrix = FeatureCounts.count_matrix
  }
}

task FastQC {
  input {
    Array[File] fastq_files_1
    Array[File] fastq_files_2
  }

  command <<<
    fastqc ~{sep=' ' fastq_files_1} ~{sep=' ' fastq_files_2} -o ./
  >>>

  output {
    Array[File] reports = glob("*.html")
  }

  runtime {
    docker: "biocontainers/fastqc:v0.11.9_cv8"
    memory: "4G"
    cpu: 2
  }
}

task STAR_Alignment {
  input {
    Array[File] fastq_files_1
    Array[File] fastq_files_2
    File reference_genome
    File gtf_annotation
    String sample_prefix
  }

  command <<<
    mkdir -p star_index
    STAR --runThreadN 4 \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles ~{reference_genome} \
         --sjdbGTFfile ~{gtf_annotation} \
         --sjdbOverhang 99

    for i in `seq 0 $((~{length(fastq_files_1)} - 1))`; do
      fq1=$(basename ~{fastq_files_1}[${i}])
      fq2=$(basename ~{fastq_files_2}[${i}])
      STAR --runThreadN 4 \
           --genomeDir star_index \
           --readFilesIn ~{fastq_files_1}[${i}] ~{fastq_files_2}[${i}] \
           --outFileNamePrefix ~{sample_prefix}_${i}_ \
           --outSAMtype BAM SortedByCoordinate
    done
  >>>

  output {
    Array[File] bam_files = glob("*.bam")
  }

  runtime {
    docker: "quay.io/biocontainers/star:2.7.9a--0"
    memory: "16G"
    cpu: 4
  }
}

task FeatureCounts {
  input {
    Array[File] bam_files
    File gtf_annotation
    String sample_prefix
  }

  command <<<
    featureCounts -T 4 -a ~{gtf_annotation} -o ~{sample_prefix}_counts.txt ~{sep=' ' bam_files}
  >>>

  output {
    File count_matrix = "~{sample_prefix}_counts.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/subread:2.0.1--h9a82719_0"
    memory: "8G"
    cpu: 4
  }
}
