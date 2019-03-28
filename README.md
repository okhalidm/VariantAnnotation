# VariantAnnotation
Annotate VCF Files

This repository will give you code to annotate VCF files using R and the ExAC API.

Start with the Challenge_data.vcf and go through the variant_annotation_script.R to give you the result of the expanded Challenge_data2.vcf file and the Annotated_VCF.csv table. The Annotated_VCF.csv table contains variant type, pertinent alternate and refrence allele frequency information, and information from the ExAC site using their REST API. 

The best way to view the Annotated_VCF.csv table would be to download it locally and open with excel or a similar program. 
A converted Annotated_VCF.xlsx has also been uploaded.

There is an additional package_install_script.R if you do not have the appropriate libraries pre-installed. Warning, runnning this script may replace packages that you already have and change your environment. You may want to install the packages independently.
