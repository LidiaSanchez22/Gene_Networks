#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

rm(list = ls())

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "C57-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

TabPhen<- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx")
rm(lnames)

#Generar sub data frame

datExprLL<- datExpr[rownames(datExpr)%in%TabPhen$ID[TabPhen$Group=="LL"],]
datExprML<- datExpr[rownames(datExpr)%in%TabPhen$ID[TabPhen$Group=="ML"],]
datExprSL<- datExpr[rownames(datExpr)%in%TabPhen$ID[TabPhen$Group=="SL"],]

#datExprLL0<- data.frame(datExprLL)
#datExprML0<- data.frame(datExprML)
#datExprSL0<- data.frame(datExprSL)



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers LL
powers = c(c(1:10), seq(from = 12, to=24, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprLL, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.0,col="red")


# Choose a set of soft-thresholding powers ML
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprML, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.0,col="red")


# Choose a set of soft-thresholding powers SL
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprSL, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.0,col="red")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


cor <- WGCNA::cor
netLL = blockwiseModules(datExprLL, power = 20,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "LLTOM", 
                       verbose = 3)

netML = blockwiseModules(datExprML, power = 30,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MLTOM", 
                       verbose = 3)

netSL = blockwiseModules(datExprSL, power = 30,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "SLTOM", 
                       verbose = 3)
cor<-stats::cor

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColorsLL = labels2colors(netLL$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netLL$dendrograms[[1]], mergedColorsLL[netLL$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColorsML = labels2colors(netML$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netML$dendrograms[[1]], mergedColorsML[netML$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColorsSL = labels2colors(netSL$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netSL$dendrograms[[1]], mergedColorsSL[netSL$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


 