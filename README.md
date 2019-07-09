# genomeseer
GenomeSeer is a collection of genomics analysis workflows that are used for detecting single nucleotide variants (SNVs), insertions/deletions (indels), copy number variants (CNVs) and translocations from RNA and DNA sequencing.  These workflows have been validated in a CLIA laboratory at UTSW

# Nextflow Workflows

## DNA Alignment

Using Module nextflow 0.31.0

```
nextflow -C nextflow.config run -w workdir alignment.nf --design design.txt --capture capture_bed_file --input dir_with_fastqs --output results_output --markdups markdup_method
```
capture_bed_file = Hyb Gene Capture Kit Bed File
markdup_methd = picard, samtools, fgbio_umi, picard_umi and none

### DNA Alignment Design File 

| SampleID | FamilyID | FqR1 | FqR2 |
|---|---|---|---|
| Sample1 | Fam1 | Sample1.R1.fastq.gz | Sample1.R2.fastq.gz |
| Sample2 | Fam1 | Sample2.R1.fastq.gz | Sample2.R2.fastq.gz |
| Sample3 | Fam2 | Sample3.R1.fastq.gz | Sample3.R2.fastq.gz |
| Sample4 | Fam2 | Sample4.R1.fastq.gz | Sample4.R2.fastq.gz |

## Tumor Only 

```
nextflow -C nextflow.config run -w workdir tumoronly.nf --design design.txt --input dir_with_bams --targetpanel capture_bed_file --output results_output
```
capture_bed_file = Hyb Gene Capture Kit Bed File

### Tumor Only Design File
| SampleID | FamilyID | BAM | GATKBAM |
|---|---|---|---|
| Sample1 | Fam1 | Sample1.dedup.bam | Sample1.final.bam |
| Sample2 | Fam1 | Sample2.dedup.bam | Sample2.final.bam |
| Sample3 | Fam2 | Sample3.dedup.bam | Sample3.final.bam |
| Sample4 | Fam2 | Sample4.dedup.bam | Sample4.final.bam |

Dedup BAMs are used for variants calling using the programs: pindel, freebayes, strelka and platypus
GATK BAMS are produced by GATK's recalibration tool and are used for variant calling using GATK.


## Somatic Variant Calling

```
nextflow -C nextflow.config run -w workdir somatic.nf --design design.txt --input dir_with_bams --output results_output
```
### Somatic Design File

| PairID | TumorID | NormalID | TumorBAM | NormalBAM | TumorGATKBAM | NormalGATKBAM |
|---|---|---|---|---|---|---|
| Subj1 | Tumor1 | Normal1 | Tumor1.dedup.bam | Normal1.dedup.bam | Tumor1.final.bam | Normal1.final.bam |
| Subj2 | Tumor2 | Normal2 | Tumor2.dedup.bam | Normal2.dedup.bam | Tumor2.final.bam | Normal2.final.bam |

Dedup BAMs are used for variants calling using the programs: pindel, freebayes, strelka, shimmer and virmid
GATK BAMS are produced by GATK's recalibration tool and are used for variant calling using MuTect2.

## RNASeq Workflow

```
nextflow -C nextflow.config run -w workdir rnaseq.nf --design design.txt --input dir_with_fastqs --output results_output --markdups markdup_method
```

markdup_methd = picard, samtools, fgbio_umi, picard_umi and none
