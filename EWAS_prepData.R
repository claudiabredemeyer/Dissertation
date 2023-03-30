
print("Loading packages...")

#suppressMessages(library("Gviz"))
#suppressMessages(library(GenomicRanges)) 
#suppressMessages(library(parallel))
#suppressMessages(library(Biobase))
#suppressMessages(library(sva))
#suppressMessages(library(limma))
#suppressMessages(library(qqman))
#suppressMessages(library(stringr))
suppressMessages(library(data.table))

args = commandArgs(T)

normalized_file <- args[1]
out_file <- args[2]

print("Loading input file...")
load(normalized_file)
print(paste("Total number of probes:",dim(betas)[1],", Total number of samples:",dim(betas)[2]))

print("Reading CrossHydridisingProbes..")
cross<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist<-cross[,1]



#[1]44210  

print("Removing CrossHydridisingProbes ...")
betas.1 <-betas[ !rownames(betas) %in% crosslist,]

print(paste(dim(betas.1)[1],"probes remained"))

print("Removing probes with SNPs...")
snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)

#[1] 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]

# [1] 10888    18
SNPlist<-snpProbes[,1]
betas.2 <- betas.1[ ! rownames(betas.1) %in% SNPlist,]

print(paste(dim(betas.2)[1],"probes remained"))

#[1] 800916    1221
#[1] 800916    1221
#[1] 800916    1221

print("Reading EPIC reference file...")
epicManifest<-fread("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7,fill=TRUE)
ANNOT<-data.frame(epicManifest)

print("Removing probes in Chr X Y...")

sexchrX<-ANNOT[which(ANNOT$CHR=="X"),]

 #[1] 19090    48
sexchrY<- ANNOT[which(ANNOT$CHR== "Y"),]
#[1] 537  48

betas.3 <- betas.2[!(rownames(betas.2) %in% sexchrX$IlmnID),]
betas.4 <- betas.3[!(rownames(betas.3) %in% sexchrY$IlmnID),]

print(paste(dim(betas.4)[1],"probes remained"))

#[1] 783656   1221
rownames(ANNOT)<-ANNOT$IlmnID

print("Removing probes with rs name...")
snps<-substr(rownames(betas.4),1,2)
rsbn<-which(snps=="rs")
if(length(rsbn)!=0)
{
  betas.5<-betas.4[-rsbn,]
}else
{
  betas.5<-betas.4
}

betas <- betas.5

print(paste(dim(betas)[1],"probes remained"))

print("Saving results...")

saveRDS(betas,file = paste0(out_file,'.prepared.rds'))

