print("Loading packages!")
suppressMessages(library("Gviz"))
suppressMessages(library(GenomicRanges))
suppressMessages(library(parallel))
suppressMessages(library(Biobase))
suppressMessages(library(sva))
suppressMessages(library(limma))
suppressMessages(library(qqman))
suppressMessages(library(stringr))
suppressMessages(library(WGCNA))
suppressMessages(library(data.table))
print("Packages are loaded.")
args<-commandArgs(TRUE)

#setwd("/mnt/data1/ClaudiaB/0603/")
#load("AdjADB.rdat")
# saveRDS(adjbetas, file="adjbetas.rds")
# 
# pheno <- read.csv("AD_CPD.csv")
#save(betas,pheno,file="ready.rdat")

print("Setting arguments...")

wd="/mnt/data1/ClaudiaB/0703"
out_pref <- ".annotated"
run_lm=1
model_lm="~Phenotype+Sex+Age+NeuNP+Sox10P+DoubleN+Plate"
num_sur=0
sva_file="ready_unadjusted_1.rdat"
#ewas_file="EWAS_RESULTS.ANNOT"

# wd <- args[1]           # Working directory
# out_pref <- args[2]     # Output files prefix
# run_lm <- args[3]       #1 means yes and 0 means No
# model_lm <- args[4]     # Model for linear regression (If run_lm set to Yes)
# num_sur <- args[5]      # Number of surrogate variables
# sva_file <- args[6]     # Input files (an .Rdata file contains 3 dataframes: 
#                                        # pheno, a dataframe contains clinical information
#                                        # betas, a matrix or dataframe contains the beta values, the variables in lm model should exist in the column names of pheno
#                                        # svs, a dataframe which its columns contain the surrugate variables)
# ewas_file <- args[7]

print("Arguments: ")
print(paste0("  Working dir= ",wd))
print(paste0("  Output files prefix= ",out_pref))
print(paste0("  Number of Surrogate variable= ",num_sur))
print(paste0("  Run linear regression to generate EWAS results= ",if(as.numeric(run_lm)==1) "Yes" else "No"))
if(as.numeric(run_lm)==1){
  print(paste0("  Name of SV file to load= ",sva_file))
}
if(as.numeric(run_lm)==1){
  print(paste0("  lm model= ",model_lm))
}
if(as.numeric(run_lm)==0){
  print(paste0("  EWAS result file to load= ",ewas_file))
}

setwd(wd)
run_lm <- as.numeric(run_lm)

if(run_lm==1){
  load(sva_file)
  
  EWAS <- function(x,model_lm,k,pheno){
    x <- as.numeric(x)
    lm_vars <- all.vars(as.formula(model_lm))
    for(j in 1:length(lm_vars)){
      
      assign(lm_vars[j], pheno[,lm_vars[j]])
    }
    model_lm <- paste0("x",model_lm)
    fit<-lm(formula = as.formula(model_lm),family = 'binomial')
    return(coef(summary(fit))[2,])
  }
  
  cl<- makeCluster(32)
  EWAS_RESULTS<-t(parApply(cl, betas , 1, EWAS, model_lm,num_sur,pheno))
  stopCluster(cl) 
  colnames(EWAS_RESULTS) <- c( "Estimate", "Std.Error"   , "t.value"  , "p.value")
  table(EWAS_RESULTS[,4]<1e-4)
  
  chisq <- qchisq(1-(EWAS_RESULTS[,4]),1)
  median(chisq)/qchisq(0.5,1)
  pdf(file = paste0(out_pref,'_qqplot_',round(median(chisq)/qchisq(0.5,1),digits = 3),'_',num_sur,'.pdf'))
  qq(EWAS_RESULTS[,4])
  dev.off()
  
  save(EWAS_RESULTS, file=paste0(out_pref,"_EWAS_RESULTS_",num_sur,".Rdata"))
}else{
  load(ewas_file)
}
p.bonf<-p.adjust(EWAS_RESULTS[,4], "bonferroni")

EWAS_RESULTS<-cbind(EWAS_RESULTS,p.bonf)
EWAS_RESULTS<-data.frame(EWAS_RESULTS)
summary(EWAS_RESULTS$p.bonf)

epicManifest<-fread("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7,stringsAsFactors = F)
ANNOT<-data.frame(epicManifest)
head(ANNOT)
x<-intersect(rownames(EWAS_RESULTS),ANNOT[,1])
ANNOT.sub<-ANNOT[which(ANNOT[,1] %in% x),]
EWAS_RESULTS<-EWAS_RESULTS[x,]

ANNOT.sub<-ANNOT.sub[order(ANNOT.sub[,1]),]
EWAS_RESULTS<-EWAS_RESULTS[order(rownames(EWAS_RESULTS)),]
identical(rownames(EWAS_RESULTS), as.character(ANNOT.sub[,1]))
#[1] TRUE

############################################################
EWAS_RESULTS.ANNOT<-cbind(EWAS_RESULTS,ANNOT.sub)
dim(EWAS_RESULTS.ANNOT)
head(EWAS_RESULTS.ANNOT)

save(EWAS_RESULTS.ANNOT, file=paste0(out_pref,"_EWAS_RESULTS_Annot_",num_sur,".Rdata"))
#EWAS_RESULTS.ANNOT$p.bh <- p.adjust(EWAS_RESULTS.ANNOT[,4], "BH")

length(which(EWAS_RESULTS.ANNOT$p.bonf<0.05))
#[1] 0
#bonf<-EWAS_RESULTS.ANNOT[which(EWAS_RESULTS.ANNOT$p.bonf<0.05),]
#table(as.character(bonf$UCSC_RefGene_Name))
#table(as.character(bonf$CHR))

############################################################
############################################################

manhattan <- data.frame("SNP"= as.character(EWAS_RESULTS.ANNOT$IlmnID),
                        "CHR"= 	as.numeric(as.character(EWAS_RESULTS.ANNOT$CHR)), 
                        "BP"= 	EWAS_RESULTS.ANNOT$MAPINFO,
                        "P"= 	EWAS_RESULTS.ANNOT[,4])
#manhattan <- manhattan[!is.na(manhattan$SNP)&!is.na(manhattan$CHR)&!is.na(manhattan$BP)&!is.na(manhattan$P),]
manhattan<-manhattan[order(manhattan$BP),]
manhattan<-manhattan[order(manhattan$CHR),]
head(manhattan)

library(qqman)

tiff(filename = paste0(out_pref,"_EWAS.MANHATTAN_",num_sur,".tiff"), units="in", width=7, height=5, res=300)
manhattan(manhattan,annotatePval=1e-5,annotateTop=FALSE,suggestiveline = FALSE, genomewideline=(-log10(1e-5)), cex.axis = 0.6 ,cex = 0.6, na.rm=na.omit,col = c("blue4", "orange3"))
dev.off()
tiff(filename = paste0(out_pref,"_EWAS.qqplot_",num_sur,".tiff"), units="in", width=3, height=3, res=300)
qq(manhattan$P, cex = 0.4, cex.axis = 0.6 , na.rm=na.omit)
dev.off()
############################################################
############################################################
EWAS_RESULTS.ANNOT$MAPINFO.1<-as.numeric(EWAS_RESULTS.ANNOT$MAPINFO)
EWAS_RESULTS.ANNOT$MAPINFO.1.1<-(EWAS_RESULTS.ANNOT$MAPINFO.1)+1

save(EWAS_RESULTS, file="EWAS_RESULTS.ANNOT")
row.names(EWAS_RESULTS.ANNOT)
save(EWAS_RESULTS.ANNOT, file="EWAS_Annotated", row.names=rownames)

IlmnID <- EWAS_RESULTS.ANNOT$IlmnID
rownames <- c(EWAS_RESULTS.ANNOT,IlmnID)

head("EWAS_Annotated")

dim(EWAS_RESULTS.ANNOT$CHR)
dim(EWAS_RESULTS.ANNOT)
#dmr <- data.frame("chrom" = 	paste0("chr", as.numeric(as.character(EWAS_RESULTS.ANNOT,|CHR))), 
dmr <- data.frame("chrom" = 	EWAS_RESULTS.ANNOT[rownames(EWAS_RESULTS.ANNOT), "CHR"],
                "start" = 	EWAS_RESULTS.ANNOT[rownames(EWAS_RESULTS.ANNOT), "MAPINFO.1"],
                  "end" 	= 	EWAS_RESULTS.ANNOT[rownames(EWAS_RESULTS.ANNOT), "MAPINFO.1.1"],
                  "pvalue" = 	EWAS_RESULTS.ANNOT[,4])

colnames(dmr) <- c("chrom", "start", "end", "pvalue")
dmr.10396 <- dmr[order(dmr[,1],dmr[,2]),]

write.table(dmr.10396, file=paste0(out_pref,"dmr_",num_sur,".txt"),sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
#comb-p pipeline -c 4 --dist 500 --seed 1.0e-3 --anno hg19 -p  /mnt/data1/Ehsan/22q11DS/dmr.22q.cnv /mnt/data1/Ehsan/22q11DS/dmr.22q.del.txt
print("All done!")