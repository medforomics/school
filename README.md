# Software for Clinical Health Omics Onocolgy Laboratories
School is a collection of genomics analysis workflows that are used for detecting single nucleotide variants (SNVs), insertions/deletions (indels), copy number variants (CNVs) and translocations from RNA and DNA sequencing.  These workflows have been validated in a CLIA laboratory at UTSW

# Nextflow Workflows

## SlideRule: DNA Workflow

using Module nextflow/20.01.0

```
nextflow -C ${nfhome}nextflow.config run -w $workdir ${nfhome}/dna.nf --input fastq --output analysis --design design.txt --seqrunid $seqrunid --pon $mutectpon --capture $capturebed --capturedir $cnvkitref --version $gittag --genome $dna_ref_path --markdups $markdup_method
```
seqrunid = Illumina Run ID, used to distrinquished repetative runs
mutectpon = capturebed = Hyb Gene Capture Kit Bed File
markdup_methd = picard, samtools, fgbio_umi, picard_umi and none

### DNA Design File 

| SampleID | TumorID | NormalID | FamilyID | FqR1 | FqR2 |
|---|---|---|---|
| Sample1 | Sample1 | Sample2 | Fam1 | Sample1.R1.fastq.gz | Sample1.R2.fastq.gz |
| Sample2 | Sample1 | Sample2 | Fam1 | Sample2.R1.fastq.gz | Sample2.R2.fastq.gz |
| Sample3 | Sample3 | Sample4 | Fam2 | Sample3.R1.fastq.gz | Sample3.R2.fastq.gz |
| Sample4 | Sample3 | Sample4 | Fam2 | Sample4.R1.fastq.gz | Sample4.R2.fastq.gz |


## Abbacus: RNASeq Workflow

```
nextflow -C nextflow.config run -w workdir rnaseq.nf --design design.txt --input dir_with_fastqs --output results_output --markdups markdup_method
```

markdup_methd = picard, samtools, fgbio_umi, picard_umi and none
