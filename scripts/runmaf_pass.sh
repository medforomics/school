#!/bin/bash

caseID="$1"

cd /archive/PHG/PHG_Clinical/cases
vepdir='/project/shared/bicf_workflow_ref/vcf2maf'
module load bedtools/2.26.0 samtools/gcc/1.8 bcftools/gcc/1.8 htslib/gcc/1.8 perl/5.28.0

vep_opt=''
tumorid=$(zgrep "#CHROM" ${caseID}/$caseID.vcf.gz|cut -f 10)
normalid=$(zgrep "#CHROM" ${caseID}/$caseID.vcf.gz|cut -f 11)
if [[ $normalid != 'NA' ]]
    then
    	vep_opt="--normal-id $normalid"
    fi
zgrep "#\|PASS" ${caseID}/$caseID.vcf.gz > ${caseID}/temp.vcf
perl ${vepdir}/vcf2maf.pl --input ${caseID}/temp.vcf --output ${caseID}/${caseID}.pass.maf --species homo_sapiens --ncbi-build GRCh38 --ref-fasta ${vepdir}/.vep/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa --filter-vcf ${vepdir}/.vep/homo_sapiens/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --cache-version 91 --vep-path ${vepdir}/variant_effect_predictor --tumor-id $tumor_id $vep_opt --custom-enst ${vepdir}/data/isoform_overrides_uniprot --custom-enst ${vepdir}/data/isoform_overrides_at_mskcc --maf-center http://www.utsouthwestern.edu/sites/genomics-molecular-pathology/ --vep-data ${vepdir}/.vep
