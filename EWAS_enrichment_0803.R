library("funr")
library(stringr)

args <- commandArgs(T)

#wd <- args[1]
input_file <- args[1]
backgr_file <- args[2]
out_pref <- args[3]

#setwd("/mnt/data1/ClaudiaB/0703/)?)

.getKEGG.1 <- function () {
  GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = "hsa", 
                                           convert = TRUE)
  isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
  GeneID.PathID <- GeneID.PathID[!isna, ]
  ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2], sep = ".")
  d <- !duplicated(ID.ID)
  GeneID.PathID <- GeneID.PathID[d, ]
  PathID.PathName <- limma::getKEGGPathwayNames(species.KEGG = "hsa", 
                                                remove.qualifier = TRUE)
  GeneID.PathID$PathwayID <- str_remove(GeneID.PathID$PathwayID,pattern = "path:")
  GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by = "PathwayID")
  kegg <- tapply(GeneID.PathID$GeneID, GeneID.PathID$PathwayID, 
                 list)
  list(idList = kegg, idTable = PathID.PathName)
}
​
gometh.1 <- function (sig.cpg, all.cpg = NULL, collection = c("GO", "KEGG"), 
                      array.type = c("450K", "EPIC"), plot.bias = FALSE, prior.prob = TRUE, 
                      anno = NULL, equiv.cpg = TRUE, fract.counts = TRUE, genomic.features = c("ALL", 
                                                                                               "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR", 
                                                                                               "ExonBnd"), sig.genes = FALSE) 
{
  array.type <- match.arg(toupper(array.type), c("450K", "EPIC"))
  collection <- match.arg(toupper(collection), c("GO", "KEGG"))
  genomic.features <- match.arg(genomic.features, c("ALL", 
                                                    "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR", 
                                                    "ExonBnd"), several.ok = TRUE)
  if (length(genomic.features) > 1 & any(grepl("ALL", genomic.features))) {
    message("All input CpGs are used for testing.")
    genomic.features <- "ALL"
  }
  if (array.type == "450K" & any(grepl("ExonBnd", genomic.features))) {
    stop("'ExonBnd' is not an annotated feature on 450K arrays,\n\n           please remove it from your genomic.feature parameter\n\n           specification.")
  }
  if (collection == "GO") {
    go <- missMethyl:::.getGO()
    result <- gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = go$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(go$idTable, result, by.x = "GOID", by.y = "row.names")
    rownames(result) <- result$GOID
  }
  else if (collection == "KEGG") {
    kegg <- .getKEGG.1()
    result <- gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = kegg$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(kegg$idTable, result, by.x = "PathwayID", 
                    by.y = "row.names")
    rownames(result) <- result$PathwayID
  }
  result[, -1]
}
####################################################################################
​
print("Loading libraries...")
suppressMessages(library("funr"))
suppressMessages(library(stringr))
suppressMessages(library(missMethyl))
suppressMessages(library("IlluminaHumanMethylationEPICmanifest"))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19 ))
suppressMessages(library(clusterProfiler))
suppressMessages(library(reshape2))
suppressMessages(library(KEGGREST))
​
args = commandArgs(T)
​
input_file <- args[1]
backgr_file <- args[2]
out_pref <- args[3]
array.type <- args[4]
print("Reading inputs...")
sigcpg <- read.csv(file = input_file,header = T,stringsAsFactors = F)[,1]
allcpg <- read.csv(file = backgr_file,header = T,stringsAsFactors = F)[,1]
print("Running Go enrichment...")
results.GO <- missMethyl::gometh(sig.cpg = sigcpg, 
                   all.cpg = allcpg, 
                   collection = "GO", 
                   array.type = array.type,
                   plot.bias = F,
                   prior.prob = T,genomic.features = "ALL",sig.genes = T)
​
results.GO <- results.GO[results.GO$N >= 10 & results.GO$N <= 2000, ]
results.GO$FDR <- p.adjust(results.GO$P.DE, method = "fdr")
​
results.GO<-results.GO[order(results.GO$P.DE),]
print("Saving GO enrichment results...")
write.csv(results.GO,file = paste0(out_pref,"_GOEnrichment.csv"))
​
​
print("Running KEGG enrichment...")
results.KEGG <- gometh.1(sig.cpg = sigcpg, 
                              all.cpg = allcpg, 
                              collection = "KEGG", 
                              array.type = array.type,
                              plot.bias = F,
                              prior.prob = T,genomic.features = "ALL",sig.genes = T)
​
results.KEGG <- results.KEGG[results.KEGG$N >= 10,]
results.KEGG$FDR <- p.adjust(results.KEGG$P.DE, method = "fdr")
​
results.KEGG<-results.KEGG[order(results.KEGG$P.DE),]
print("Saving KEGG enrichment results...")
write.csv(results.KEGG,file = paste0(out_pref,"_KeggEnrichment.csv"))


print("Reading inputs...")
#source(paste0(str_remove(sys.script(),pattern="EWAS_enrichment.R"),"crystalmethEPIC.R"))
source("crystalmethEPIC.R")

sigcpg <- read.csv(file = input_file,header = T,stringsAsFactors = F)[,1]
allcpg <- read.csv(file = backgr_file,header = T,stringsAsFactors = F)[,1]


print("KEGG enrichment...")
kegg <- crystalmeth(sig.cpg = sigcpg,
                    all.cpg = allcpg,
                    collection = "KEGG", array.type = "EPIC" )

print("Saving KEGG results...")                    
kegg<-kegg[order(go$P.DE),]
write.csv(kegg,file = paste0(out_pref,"_KeggEnrichment.csv"))

print("All Done!")

print("GO enrichment...")
go <- crystalmeth(sig.cpg = sigcpg,
                  all.cpg = allcpg,
                  collection = "GO",array.type = "EPIC" )

print("Saving GO results...")
go<-go[order(go$P.DE),]
#head(go[,1:6], 20)



print("All Done!")
