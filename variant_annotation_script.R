
#Set your working directory to when you have saved the "Challenge_data.vcf" file
#Call in library for Variant Annotation
#Read in the variant data and expand out multiallelic sites
#Write the expanded set as "Challenge_data2.vcf" which has been provided in this github

library(VariantAnnotation)
vcf1 <- readVcf("Challenge_data.vcf", "hg19")
vcf2 <- expand(vcf1)
writeVcf(vcf2,"Challenge_data2.vcf")

#Read in the expanded VCF file

vcf1 <- readVcf("Challenge_data2.vcf", "hg19")

#obtain the header and genotype information of the vcf file to help in downstream analysis. The two lines below have been commented out but may be useful to those who which to look over that information.

#head(vcf1)
#geno(vcf1)


#Identify the different variant types in the expanded VCF file
Deletions <- isDeletion(vcf1)
Indels <- isIndel(vcf1)
Insertions <- isInsertion(vcf1)
Substitutions <- isSubstitution(vcf1)
Delins <- isDelins(vcf1)
SNVs <- isSNV(vcf1)
Transitions<- isTransition(vcf1)

#obtain the vcf header information and corresponding data of the VCF file
info_vcf1 <- vcf1@info

#Column bind the vcf data
info_vcf1.1 <- cbind(info_vcf1, Deletions, Indels, Insertions, Substitutions, Delins, SNVs, Transitions)

#Call in vcfR package to manipulate the VCF datafram
library(vcfR)

#obtain the fixed and genotype data for the vcf - this will allow you to work with the different samples in the vcf file

y <- read.vcfR("Challenge_data2.vcf")
y@fix[,c(1:5)]
q <- as.data.frame(y@fix[,c(1:5)])
q1 <-as.data.frame(y@gt)

#Call in the splitstackshape library and expand out the normal and Vaf5 data and give appropriate headers

library(splitstackshape)
q1 <-cbind(q,q1)

q1.1 <- cSplit(q1,c("normal","vaf5"), ":")
colnames(q1.1) <- c("CHROM","POS","ID","REF","ALT", "FORMAT", "normal_GT", "normal_GQ", "normal_DP", "normal_DPR", "normal_RO","normal_QR", "normal_A0","normal_QA","vaf5_GT", "vaf5_GQ", "vaf5_DP", "vaf5_DPR", "vaf5_RO","vaf5_QR", "vaf5_A0","vaf5_QA")

#Calculate frequencies of the reference allele (RO) and the alternate allele (A0) against the read depth (DP) for the different samples in the vcf data

q1.1$normal_ref_freq <- as.numeric(as.character(q1.1$normal_RO))/as.numeric(as.character(q1.1$normal_DP))
q1.1$vaf5_ref_freq <- as.numeric(as.character(q1.1$vaf5_RO))/as.numeric(as.character(q1.1$vaf5_DP))

q1.1$normal_alt_freq <- as.numeric(as.character(q1.1$normal_A0))/as.numeric(as.character(q1.1$normal_DP))
q1.1$vaf5_alt_freq <- as.numeric(as.character(q1.1$vaf5_A0))/as.numeric(as.character(q1.1$vaf5_DP))


#Pull ExAC data by creating a list of hyperlinks to ping against ExAC and pull the corresponding json files

z <- paste("exac.hms.harvard.edu/rest/variant/variant/",q$CHROM,"-",q$POS,"-",q$REF,"-",q$ALT,sep="")

#Call in libraries to pull the data from the ExAC site using the links created above

library(lubridate)
library(httr)
library(jsonlite)

options(stringsAsFactors = FALSE)

#Create a for loop to get the allele frequencies
t <- NULL
x <- seq(1,length(z),1)
for(i in x){url <- z[i]
raw.result <- GET(url = url)
this.raw.content <- rawToChar(raw.result$content)
this.content <- fromJSON(this.raw.content)
s<-this.content$allele_freq
if(length(s)==0) {
  s = "NA"
} else {
  s<-this.content$allele_freq
}
t <- rbind(t,s)
}

t1 <- as.data.frame(t)
colnames(t1) <- "allele_freq"

#Create a for loop to pull the rsid

u <- NULL
x <- seq(1,length(z),1)
for(i in x){url <- z[i]
raw.result <- GET(url = url)
this.raw.content <- rawToChar(raw.result$content)
this.content <- fromJSON(this.raw.content)
s<-this.content$rsid
if(length(s)==0) {
  s = "NA"
} else {
  s<-this.content$rsid
}
u <- rbind(u,s)
}

u1 <- as.data.frame(u)
colnames(u1) <- "rsid"

#Create a for loop to pull the gene symbol
v <- NULL
x <- seq(1,length(z),1)
for(i in x){url <- z[i]
raw.result <- GET(url = url)
this.raw.content <- rawToChar(raw.result$content)
this.content <- fromJSON(this.raw.content)
s<-unique(this.content$vep_annotations$SYMBOL)
if(length(s)==0) {
  s = "NA"
} else {
  s<-unique(this.content$vep_annotations$SYMBOL)
}
v <- rbind(v,s)
}

v1 <- as.data.frame(v)
colnames(v1) <- "Gene_Symbol"

#Create a for loop to pull the major_consequence

w <- NULL
x <- seq(1,length(z),1)
for(i in x){url <- z[i]
raw.result <- GET(url = url)
this.raw.content <- rawToChar(raw.result$content)
this.content <- fromJSON(this.raw.content)
s<-unique(this.content$vep_annotations$major_consequence)
if(length(s)==0) {
  s = "NA"
} else {
  s<-unique(this.content$vep_annotations$major_consequence)
}
w <- rbind(w,s)
}

w1 <- as.data.frame(w)
colnames(w1) <- "Major_Consequence"

#Create a for loop to pull the existing variations
n <- NULL
x <- seq(1,length(z),1)
for(i in x){url <- z[i]
raw.result <- GET(url = url)
this.raw.content <- rawToChar(raw.result$content)
this.content <- fromJSON(this.raw.content)
s<-unique(this.content$vep_annotations$Existing_variation)
if(length(s)==0) {
  s = "NA"
} else {
  s<-unique(this.content$vep_annotations$Existing_variation)
}
n <- rbind(n,s)
}

n1 <- as.data.frame(n)
colnames(n1) <- "Existing_variation"

#Combine the data frames together and extract the corresponding data for the variant, variant type, and frequencies
t2 <- cbind(info_vcf1.1, q1.1, t1, u1,v1,w1,n1)
t3 <- t2[c(51:55,44:50,33,59,61,63,67,69,71,73:81)]

#write a csv table of the annotated variants which is uploaded onto this github repository
write.csv(t3,"Annotated_VCF.csv")
