from __future__ import print_function
from pysam import VariantFile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--caseID', help='case ID for input vcf')
args = parser.parse_args()

infile=args.caseID + ".vcf.gz"
outfile=args.caseID + ".filtered.vcf"
errfile=args.caseID + ".fusions.vcf"
vcf_in = VariantFile(infile)  # auto-detect input format
vcf_out = VariantFile(outfile, 'w', header=vcf_in.header)
vcf_err = VariantFile(errfile, 'w', header=vcf_in.header)
vcf_out.close()
vcf_err.close()

with open(outfile, 'w') as out, open(errfile, 'w') as err:
	for rec in vcf_in.fetch():
		record=str(rec).rstrip()
		line=record.split("\t")
		infoArray=line[7].split(";")

		effects=0
		fusions=0
		for field in infoArray:
			# check to see if the section starts with "ANN:"
			if field[0:4] == "ANN=":
				annotation=field[4:]
				annArray=annotation.split(",")
				for SNPeff in annArray:
					SNPeffArray=SNPeff.split("|")
					if SNPeffArray[2] == "HIGH" or SNPeffArray[2] == "MODERATE":
						effects += 1
						if "gene_fusion" in SNPeffArray[1]:
							fusions += 1
		outline='\t'.join(line)	
		if effects == fusions and fusions > 0:
			print(outline, file=err)
		else:
			print(outline, file=out)

vcf_in.close()
out.close()
err.close()
