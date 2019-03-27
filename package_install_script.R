#Package installations script if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation", version = "3.8")

install.packages(c("vcfR","lubridate","httr","jsonlite","splitstackshape")
