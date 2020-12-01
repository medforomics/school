task trimseq {
  input {
    File fq1
    File fq2
    String sampleid
  }
  command {
    bash ${repoDir}/process_scripts/preproc_fastq/trimgalore.sh -f -p ${sampleid} ${fq1} ${fq2}
    
  }
  output {
    File tfq1 = "${sampleid}.trim.R1.fastq.gz"
    File tfq2 = "${sampleid}.trim.R2.fastq.gz"
    File treport = "${sampleid}.trimreport.txt"
  }
  runtime {
    docker: "goalconsortium/trim_galore:1.0.7"
    cpu: 1
    memory: "4 GB"
  }
}

task dna_align {
  input {
    File fq1
    File fq2
    
  }

}

workflow dna {
  inputs {
    File FqR1
    File FqR2
    String CID
    String SID
    String paneldir
    String refdir
    
  }
  call trimseq {
    input:
      fq1=FqR1
      fq2=FqR1
      sampleid=SID
  }
}
