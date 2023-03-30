
print("Loading packages...")
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))

args<-commandArgs(TRUE)
print("Reading inputs...")
beta_file <- args[1] # RDS file
pheno_file <- args[2] # csv file
variables_numeric <- args[3]
variables_numeric <- str_split(variables_numeric,pattern = ',',simplify = T)[1,]
variables_factor <- args[4]
variables_factor <- str_split(variables_factor,pattern = ',',simplify = T)[1,]
out_prefix <- args[5]

betas <- readRDS(beta_file)
pheno <- read.csv(pheno_file,row.names =NULL,stringsAsFactors = F)
pheno <- pheno[,c(variables_numeric, variables_factor)]
#if("sentrixid" %in% tolower(colnames(pheno))){
#  index <- which(tolower(colnames(pheno)) %in% "sentrixid")
#  pheno[,index] <- str_split(rownames(pheno),pattern = "_",simplify = T)[,1]
#}
print("Calculating PCs...")
pca <- prcomp(t(betas),center = T,scale. = T)
pca1 <- pca$x[,1:10]
pca1 <- data.frame(pca1)

cor_ <- matrix(data = NA,nrow =10,ncol = length(c(variables_numeric,variables_factor)) )
colnames(cor_) <- c(variables_numeric,variables_factor)
rownames(cor_) <- colnames(pca1)[1:10]

cor_pval <- cor_
print("Calculating correlations...")
for (i in 1:10) {
  for (j in 1:length(c(variables_numeric,variables_factor))) {
    var_ <- c(variables_numeric,variables_factor)[j]
    
    if(var_ %in% variables_factor){
      res1<-cor.test(as.numeric(pca1[,i]),as.numeric(as.factor(pheno[,var_])), method="spearman",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }else{
      res1<-cor.test(as.numeric(pca1[,i]),as.numeric(pheno[,var_]), method="pearson",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }
  }
}

textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix) = dim(cor_)
print("Saving Plot...")
pdf(file = paste0(out_prefix,".PCA.pdf"))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = rownames(cor_),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("PCA Analysis"))

dev.off()

