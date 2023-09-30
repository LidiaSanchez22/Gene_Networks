#Cargar las libraries necesarias
rm(list = ls())
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(xlsx)
library(readxl)
library(dplyr)
library(ggvenn)
library(ggVennDiagram)


#
# ===== Code chunk 1 Load Data =====
#


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/OBESITAT";
setwd(workingDir); 

# Load aligment data
load(file = "C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/RSTUDIO-EXPORT/aligned_JZJOS_02.Rdata" )

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
COLNAMES <- c(substr(colnames(aligned$counts),1,nchar(colnames(aligned$counts))-4))
DF_counts <- data.frame(aligned$counts) %>% `colnames<-`(COLNAMES)
IlluminaID <- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/RSTUDIO-EXPORT/JIMENEZJOS_02.xls",sheet = "IDs_WAT")
CNAGtoLAB_ID <- data.frame(IlluminaID$...10, IlluminaID$...14, IlluminaID$...9, IlluminaID$...2,  IlluminaID$...12,  IlluminaID$...13) %>% 
  .[c(-1,-2),] %>% `colnames<-`(c("CNAG_ID", "LAB_ID", "Machine", "Lane", "Group", "Generation")) %>% na.omit()
rm(IlluminaID)
data <- aligned$counts 
for ( col in 1:ncol(data)){
  colnames(data)[col] <-  sub(".bam", "", colnames(data)[col])
}
data <- data.frame(data)
data <- data[,(colnames(data) %in% CNAGtoLAB_ID$CNAG_ID)]

#
# ===== Code chunk 2 Quality check and  outlier detection =====
#


gsg <- goodSamplesGenes(t(data)) # input is a matrix or data frame in which columns are genes and rows ar samples
summary(gsg)
# to detect whether there are outliers in our data
gsg$allOK  # if TRUE -> all the samples and the genes passed, so none of them are outliers or missing entries.
# if FLASE -> either the genes or samples are outliers or massing value. we need to find out and exclude them

table(gsg$goodGenes) # -> In this case 1114 genes are FALSE and have to be removed
#FALSE  TRUE 
#1114 26065
table(gsg$goodSamples) # In this case all samples are TRUE, so passed the test.
#TRUE 
#316

#
# ===== Code chunk 3 - Remove genes which are detected as outliers =====
#

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) {
    printFlush(paste("Removing genes:", paste(names(data)[!gsg$goodGenes], collapse = ", ")));
  }
     
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", ")));
}
# Keep only those genes and samples that are labeled as GOOD
data <- data[gsg$goodGenes == TRUE,]
data <- data[gsg$goodSamples == TRUE,]

#
# ===== Code chunk 4 - Detect outlier samples by hierarchical clustering =====
# I PREFER TO SKIP THIS AND USE PCA SINCE I CAN HAVE A LOOK AT BATCH EFFECTS ALSO

# Detectio outlier samples - hierarchical clustering - method 1
sampleTree = hclust(dist(t(data)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 5)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
data = data[keepSamples, ]
nGenes = ncol(data)
nSamples = nrow(data)

#
# =====  Code chunk 5 - Detect outlier samples by PCA =====
#

### NOTE: If there are batch effects observed, correct for them before moving ahead
# I'll apply PCA to see whether there are batch effects due to machine 
MACHINE <- CNAGtoLAB_ID$Machine[CNAGtoLAB_ID$CNAG_ID %in% colnames(data)]
LANE <- CNAGtoLAB_ID$Lane[CNAGtoLAB_ID$CNAG_ID %in% colnames(data)]

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(x=PC1, y=PC2, color=MACHINE, shape = LANE)) +
  geom_point() +
  #geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# from the graph we can see that there are some clustering 
SampltoRem <- c(rownames(pca.dat[pca.dat$PC1>5000000,]), rownames(pca.dat[pca.dat$PC2<-1000000,]))
data <- data[,!(colnames(data) %in% SampltoRem)]
rm(pca, pca.dat, pca.var,MACHINE, LANE)
MACHINE <- CNAGtoLAB_ID$Machine[CNAGtoLAB_ID$CNAG_ID %in% colnames(data)]
LANE <- CNAGtoLAB_ID$Lane[CNAGtoLAB_ID$CNAG_ID %in% colnames(data)]


pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2, color=MACHINE, shape = LANE)) +
  geom_point() +
  #geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
  ggtitle("PCA for batch effects")

rm(pca, pca.dat, pca.var,MACHINE, LANE)
#Sin batch effect en lane ni por machine (Netejar el gràfic sense noms) - Justifica no correccions per batch effect

#
# ===== Code chunk 6 - Normalization  =====
#
  
# WGCNA requires data to be normalized. Since we have raw counts, we need to normalize them 
# the creator the packages recommended to use variance stabilizing transformation from the DESeq2 package.

# create a deseq2 dataset
phenoData <- CNAGtoLAB_ID[CNAGtoLAB_ID$CNAG_ID %in% colnames(data),c(1,2,5,6)]
  
# making the rownames and column names identical
phenoData<-phenoData[order(match(phenoData[[1]], colnames(data))), ] %>%
    remove_rownames() %>%
    column_to_rownames(var = 'CNAG_ID')
  
all(rownames(phenoData) %in% colnames(data)) # has to be TRUE 
all(rownames(phenoData) == colnames(data))  # has to be TRUE 
  
# create dds dataset
  dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = phenoData,
                                design = ~ 1) # not specifying model because we need DESeq to perform variance stabilizing transformation.
  ddsColl <- collapseReplicates(dds, dds$LAB_ID)

## remove all genes with counts < 5 in more than 75% of samples (88*0.75=66)
## suggested by WGCNA on RNAseq FAQ
dds75 <- ddsColl[rowSums(counts(ddsColl) >= 5) >= 66,] #son estos numeros como podrian ser otros, es un poco alternativo
  nrow(dds75) # genes
  
# perform variance stabilization
  dds_norm <- vst(dds75)
  
# get normalized counts
  norm.counts <- assay(dds_norm) %>% 
    t() # transpose is needed since the next function requires the rows to be the samples and the columns to be the gene IDs
 
  
  #
  # ===== Code chunk 7 - Traits generation =====
  #
  
  TabPhen<- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx")
  
  #quitar duplicados del phenodata (88columnes)
  #unique(phenoData, incomparables = FALSE, fromLast = FALSE,nmax = NA)  #buscar funciones para data frame
  
  #select rows with unique values based on LAB_ID column only
  phenoDataO <- phenoData %>% distinct(LAB_ID, .keep_all=TRUE)
  
  #phenoData = TabPhen %>% dplyr::select(ID, Group, Generation) %>% distinct( )
  
  #phenoData <- TabPhen[!duplicated(TabPhen), ]
  
  #Reordenar las columnas para que la 1a sea el Lab_ID dentro del phenoDataO
  phenoDataX <- phenoDataO[,-1]
  rownames(phenoDataX) <- phenoDataO[,1]
  
  all(rownames(norm.counts) %in% rownames(phenoDataX)) # has to be TRUE
  
  #Ordenar las filas para que esten en el mismo orden en phenoDataO y norm.counts
  
  phenoDataX <- phenoDataX[rownames(norm.counts),]
  
  all(rownames(norm.counts) == rownames(phenoDataX))  # has to be TRUE

  #datTraits <- data.frame(cbind(rownames(phenoDataX),phenoDataX$Generation)) %>%
   # `colnames<-`(c("LAB_ID", "Generation"))
 # row.names(datTraits) <- datTraits[,1]
  #datTraits[,1] <- NULL
  
 # for (irow in c(1:length(rownames(phenoDataX)))) {
  #  if (rownames(phenoDataX)[irow] %in% TabPhen$ID){datTraits[irow,c(2:24)]<-TabPhen[which(TabPhen$ID == rownames(phenoDataX)[irow]), c(2,14:35) ]}
 # }
  
  #library("")#buscar library
  #datTraits<- datTraits[complete.cases(datTraits[ , 5:11]),]
  datTraits<- TabPhen[, c(1:3, 14:27, 32:34)]

  
  datExpr <- norm.counts[rownames(norm.counts) %in% datTraits$ID,]
  collectGarbage();
  
  rm(norm.counts, data, dds, dds_norm, dds75, gsg, CNAGtoLAB_ID, DF_counts)
  
  #Convertimos todas las variables dentro de datTraits a númericas
  
  #datTraits$LAB_ID<-as.factor(datTraits$LAB_ID)
  datTraits$Group<-as.factor(datTraits$Group)
  datTraits$Generation<-as.factor(datTraits$Generation)
  #datTraits$`BW.SAC`<-as.numeric(datTraits$`BW.SAC`)
  datTraits$`GLU.SAC`<-as.numeric(datTraits$`GLU.SAC`)
  datTraits$`INS.SAC`<-as.numeric(datTraits$`INS.SAC`)
  datTraits$`INS0MIN`<-as.numeric(datTraits$`INS0MIN`)
  datTraits$`INS30MIN`<-as.numeric(datTraits$`INS30MIN`)
  datTraits$`GLU.0`<-as.numeric(datTraits$`GLU.0`)
  datTraits$`GLU.30`<-as.numeric(datTraits$`GLU.30`)
  datTraits$`RATIO INS/GLU0`<-as.numeric(datTraits$`RATIO INS/GLU0`)
  datTraits$`RATIO INS/GLU 30`<-as.numeric(datTraits$`RATIO INS/GLU 30`)
  datTraits$`D.INS`<-as.numeric(datTraits$`D.INS`)
  datTraits$`D.GLU`<-as.numeric(datTraits$`D.GLU`)
  datTraits$LIVER<-as.numeric(datTraits$LIVER)
  #datTraits$Ewat<-as.numeric(datTraits$Ewat)
  #datTraits$Iwat<-as.numeric(datTraits$Iwat)
  datTraits$`TAG.SAC`<-as.numeric(datTraits$`TAG.SAC`)
  datTraits$HOMA<-as.numeric(datTraits$HOMA.IR)
  datTraits$AUC<-as.numeric(datTraits$AUC)
  #datTraits$`GR 0-7`<-as.numeric(datTraits$`GR 0-7`)
  datTraits$`percentage.liver`<-as.numeric(datTraits$`percentage.liver`)
  #datTraits$`% eWAT`<-as.numeric(datTraits$`% eWAT`)
  #datTraits$`% iWAT`<-as.numeric(datTraits$`% iWAT`)
  #datTraits$...35<-as.numeric(datTraits$...35)
  
  
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(as.matrix(datTraits[,c(1:3,14:24)]), signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits[,c(1:3,14:24)]),
                      main = "Sample dendrogram and trait heatmap")
 
  #el error se debe a que datExpr solo tiene 19 columnas, por lo que en dataTraits, tenemos 26, para obtener el hetmap, solo podemos seleccionar hasta la 19
  
  
  
  #
  # ===== Code chunk 8 - Save =====
  #
  
  save(datExpr, datTraits, file = "C57-01-dataInput.RData")
 