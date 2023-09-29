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
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/Documents/JIMENEZJOS_02/Res_WGCNA_auto";
setwd(workingDir); 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "C57-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

# Load network data saved in the second part.
lnamesLL = load(file = "MouseLL-networkConstruction-auto.RData");
lnamesML = load(file = "MouseML-networkConstruction-auto.RData");
lnamesSL = load(file = "MouseSL-networkConstruction-auto.RData");
rm(lnames, lnamesLL, lnamesML, lnamesSL)

#Generar sub data frame

TabPhen<- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx")

datExprLL<- datExpr[rownames(datExpr)%in%TabPhen$ID[TabPhen$Group=="LL"],]
datExprML<- datExpr[rownames(datExpr)%in%TabPhen$ID[TabPhen$Group=="ML"],]
datExprSL<- datExpr[rownames(datExpr)%in%TabPhen$ID[TabPhen$Group=="SL"],]

nGenesLL = ncol(datExprLL)
nSamplesLL = nrow(datExprLL)

nGenesML = ncol(datExprML)
nSamplesML = nrow(datExprML)

nGenesSL = ncol(datExprSL)
nSamplesSL = nrow(datExprSL)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOMLL = 1-TOMsimilarityFromExpr(datExprLL, power = 20);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOMLL = dissTOMLL^21;
# Set diagonal to NA for a nicer plot
diag(plotTOMLL) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOMLL, geneTreeLL, moduleColorsLL, main = "Network heatmap plot LL, all genes")


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOMML = 1-TOMsimilarityFromExpr(datExprML, power = 30);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOMML = dissTOMML^31;
# Set diagonal to NA for a nicer plot
diag(plotTOMML) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOMML, geneTreeML, moduleColorsML, main = "Network heatmap plot ML, all genes")



# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOMSL = 1-TOMsimilarityFromExpr(datExprSL, power = 30);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOMSL = dissTOMSL^31;
# Set diagonal to NA for a nicer plot
diag(plotTOMSL) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOMSL, geneTreeSL, moduleColorsSL, main = "Network heatmap plot SL, all genes")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
selectLL = sample(nGenesLL, size = nSelect);
selectTOMLL = dissTOMLL[selectLL, selectLL];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTreeLL = hclust(as.dist(selectTOMLL), method = "average")
selectColorsLL = moduleColorsLL[selectLL];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDissLL = selectTOMLL^21;
diag(plotDissLL) = NA;
TOMplot(plotDissLL, selectTreeLL, selectColorsLL, main = "Network heatmap plot, selected genes LL")


nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
selectML = sample(nGenesML, size = nSelect);
selectTOMML = dissTOMML[selectML, selectML];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTreeML = hclust(as.dist(selectTOMML), method = "average")
selectColorsML = moduleColorsML[selectML];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDissML = selectTOMML^7;
diag(plotDissML) = NA;
TOMplot(plotDissML, selectTreeML, selectColorsML, main = "Network heatmap plot, selected genes ML")


nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
selectSL = sample(nGenesSL, size = nSelect);
selectTOMSL = dissTOMSL[selectSL, selectSL];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTreeSL = hclust(as.dist(selectTOMSL), method = "average")
selectColorsSL = moduleColorsSL[selectSL];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDissSL = selectTOMSL^7;
diag(plotDissSL) = NA;
TOMplot(plotDissSL, selectTreeSL, selectColorsSL, main = "Network heatmap plot, selected genes SL")


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Recalculate module eigengenes
MEsLL = moduleEigengenes(datExprLL, moduleColorsLL)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(subdatTraitsLL$percentage.WAT);
names(weight) = "% WAT"
# Add the weight to existing module eigengenes
METLL = orderMEs(cbind(MEsLL, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(METLL, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Recalculate module eigengenes
MEsML = moduleEigengenes(datExprML, moduleColorsML)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(subdatTraitsML$percentage.WAT);
names(weight) = "% WAT"
# Add the weight to existing module eigengenes
METML = orderMEs(cbind(MEsML, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(METML, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


# Recalculate module eigengenes
MEsSL = moduleEigengenes(datExprSL, moduleColorsSL)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(subdatTraitsSL$percentage.WAT);
names(weight) = "% WAT"
# Add the weight to existing module eigengenes
METSL = orderMEs(cbind(MEsSL, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(METSL, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(METLL, "Eigengene dendrogram LL", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(METLL, "Eigengene adjacency heatmap LL", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(METML, "Eigengene dendrogramML", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(METML, "Eigengene adjacency heatmapML", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(METSL, "Eigengene dendrogramSL", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(METSL, "Eigengene adjacency heatmap SL", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
