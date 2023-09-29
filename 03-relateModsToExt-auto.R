#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/OBESITAT";
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
lnamesLL

lnamesML = load(file = "MouseML-networkConstruction-auto.RData");
lnamesML

lnamesSL = load(file = "MouseSL-networkConstruction-auto.RData");
lnamesSL
rm(lnames, lnamesLL, lnamesML, lnamesSL)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenesLL = ncol(datExprLL);
nSamplesLL = nrow(datExprLL);

nGenesML = ncol(datExprML);
nSamplesML = nrow(datExprML);

nGenesSL = ncol(datExprSL);
nSamplesSL = nrow(datExprSL);

subdatTraitsLL<- datTraits[datTraits$Group=="LL",]
subdatTraitsML<- datTraits[datTraits$Group=="ML",]
subdatTraitsSL<- datTraits[datTraits$Group=="SL",]

# Recalculate MEs with color labels
MEs0LL = moduleEigengenes(datExprLL, moduleColorsLL)$eigengenes
MEsLL = orderMEs(MEs0LL)
subdatTraitsLL$`GLU.30` <- as.numeric(subdatTraitsLL$`GLU.30`)
subdatTraitsLL$`INS.SAC` <- as.numeric(subdatTraitsLL$`INS.SAC`)
subdatTraitsLL <- subdatTraitsLL[subdatTraitsLL$ID%in%rownames(MEsLL),c(4:ncol(subdatTraitsLL))] # extract only traits/clinical variables
moduleTraitCorLL = cor(MEsLL, subdatTraitsLL, use = "p"); # pearson correlation
moduleTraitPvalueLL = corPvalueStudent(moduleTraitCorLL, nSamplesLL);


MEs0ML = moduleEigengenes(datExprML, moduleColorsML)$eigengenes
MEsML = orderMEs(MEs0ML)
subdatTraitsML$`GLU.30` <- as.numeric(subdatTraitsML$`GLU.30`)
subdatTraitsML$`INS.SAC` <- as.numeric(subdatTraitsML$`INS.SAC`)
subdatTraitsML <- subdatTraitsML[subdatTraitsML$ID%in%rownames(MEsML),c(4:ncol(subdatTraitsML))] # extract only traits/clinical variables
moduleTraitCorML = cor(MEsML, subdatTraitsML, use = "p"); # pearson correlation
moduleTraitPvalueML = corPvalueStudent(moduleTraitCorML, nSamplesML);

MEs0SL = moduleEigengenes(datExprSL, moduleColorsSL)$eigengenes
MEsSL = orderMEs(MEs0SL)
subdatTraitsSL$`GLU.30` <- as.numeric(subdatTraitsSL$`GLU.30`)
subdatTraitsSL$`INS.SAC` <- as.numeric(subdatTraitsSL$`INS.SAC`)
subdatTraitsSL <- subdatTraitsSL[subdatTraitsSL$ID%in%rownames(MEsSL),c(4:ncol(subdatTraitsSL))] # extract only traits/clinical variables
moduleTraitCorSL = cor(MEsSL, subdatTraitsSL, use = "p"); # pearson correlation
moduleTraitPvalueSL = corPvalueStudent(moduleTraitCorSL, nSamplesSL);

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCorLL, 2), "\n(",
                    signif(moduleTraitPvalueLL, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorLL)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCorLL,
               xLabels = names(subdatTraitsLL),
               yLabels = names(MEsLL),
               ySymbols = names(MEsLL),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships LL"))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCorML, 2), "\n(",
                    signif(moduleTraitPvalueML, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorML)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCorML,
               xLabels = names(subdatTraitsML),
               yLabels = names(MEsML),
               ySymbols = names(MEsML),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships ML"))


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCorSL, 2), "\n(",
                    signif(moduleTraitPvalueSL, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorSL)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCorSL,
               xLabels = names(subdatTraitsSL),
               yLabels = names(MEsSL),
               ySymbols = names(MEsSL),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships SL"))

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

##Large Lite LL##
# Define variable weight containing the weight column of datTrait
percentageWATLL= as.data.frame(subdatTraitsLL$percentage.WAT);
names(percentageWATLL) = "% WAT LL"
# names (colors) of the modules
modNamesLL = substring(names(MEsLL), 3)

geneModuleMembershipLL = as.data.frame(cor(datExprLL, MEsLL, use = "p"));
MMPvalueLL = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipLL), nSamplesLL));

names(geneModuleMembershipLL) = paste("MM", modNamesLL, sep="");
names(MMPvalueLL) = paste("p.MM", modNamesLL, sep="");

geneTraitSignificanceLL = as.data.frame(cor(datExprLL, percentageWATLL, use = "p"));
GSPvalueLL = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceLL), nSamplesLL));

names(geneTraitSignificanceLL) = paste("GS.", names(percentageWATLL), sep="");
names(GSPvalueLL) = paste("p.GS.", names(percentageWATLL), sep="");


##Medium Lite ML##
# Define variable weight containing the weight column of datTrait
percentageWATML= as.data.frame(subdatTraitsML$percentage.WAT);
names(percentageWATML) = "% WAT ML"
# names (colors) of the modules
modNamesML = substring(names(MEsML), 3)

geneModuleMembershipML = as.data.frame(cor(datExprML, MEsML, use = "p"));
MMPvalueML = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipML), nSamplesML));

names(geneModuleMembershipML) = paste("MM", modNamesML, sep="");
names(MMPvalueML) = paste("p.MM", modNamesML, sep="");

geneTraitSignificanceML = as.data.frame(cor(datExprML, percentageWATML, use = "p"));
GSPvalueML = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceML), nSamplesML));

names(geneTraitSignificanceML) = paste("GS.", names(percentageWATML), sep="");
names(GSPvalueML) = paste("p.GS.", names(percentageWATML), sep="");


##Medium Lite SL##
# Define variable weight containing the weight column of datTrait
percentageWATSL= as.data.frame(subdatTraitsSL$percentage.WAT);
names(percentageWATSL) = "% WAT SL"
# names (colors) of the modules
modNamesSL = substring(names(MEsSL), 3)

geneModuleMembershipSL = as.data.frame(cor(datExprSL, MEsSL, use = "p"));
MMPvalueSL = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipSL), nSamplesSL));

names(geneModuleMembershipSL) = paste("MM", modNamesSL, sep="");
names(MMPvalueSL) = paste("p.MM", modNamesSL, sep="");

geneTraitSignificanceSL = as.data.frame(cor(datExprSL, percentageWATSL, use = "p"));
GSPvalueSL = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceSL), nSamplesSL));

names(geneTraitSignificanceSL) = paste("GS.", names(percentageWATSL), sep="");
names(GSPvalueSL) = paste("p.GS.", names(percentageWATSL), sep="");

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#Correlation cloud for LL with most significant modules#

moduleLL = "darkmagenta"
column = match(moduleLL, modNamesLL);
moduleGenesLL = moduleColorsLL==moduleLL;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipLL[moduleGenesLL, column]),
                   abs(geneTraitSignificanceLL[moduleGenesLL, 1]),
                   xlab = paste("Module Membership in", moduleLL, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleLL)

moduleLL = "pink"
column = match(moduleLL, modNamesLL);
moduleGenesLL = moduleColorsLL==moduleLL;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipLL[moduleGenesLL, column]),
                   abs(geneTraitSignificanceLL[moduleGenesLL, 1]),
                   xlab = paste("Module Membership in", moduleLL, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleLL)

moduleLL = "darkslateblue"
column = match(moduleLL, modNamesLL);
moduleGenesLL = moduleColorsLL==moduleLL;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipLL[moduleGenesLL, column]),
                   abs(geneTraitSignificanceLL[moduleGenesLL, 1]),
                   xlab = paste("Module Membership in", moduleLL, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleLL)



#Correlation cloud for ML with most significant modules#

moduleML = "pink"
column = match(moduleML, modNamesML);
moduleGenesML = moduleColorsML==moduleML;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipML[moduleGenesML, column]),
                   abs(geneTraitSignificanceML[moduleGenesML, 1]),
                   xlab = paste("Module Membership in", moduleML, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleML)

moduleML = "grey60"
column = match(moduleML, modNamesML);
moduleGenesML = moduleColorsML==moduleML;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipML[moduleGenesML, column]),
                   abs(geneTraitSignificanceML[moduleGenesML, 1]),
                   xlab = paste("Module Membership in", moduleML, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleML)


moduleML = "black"
column = match(moduleML, modNamesML);
moduleGenesML = moduleColorsML==moduleML;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipML[moduleGenesML, column]),
                   abs(geneTraitSignificanceML[moduleGenesML, 1]),
                   xlab = paste("Module Membership in", moduleML, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleML)


moduleML = "midnightblue"
column = match(moduleML, modNamesML);
moduleGenesML = moduleColorsML==moduleML;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipML[moduleGenesML, column]),
                   abs(geneTraitSignificanceML[moduleGenesML, 1]),
                   xlab = paste("Module Membership in", moduleML, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleML)


#Correlation cloud for SL with most significant modules#

moduleSL = "turquoise"
column = match(moduleSL, modNamesSL);
moduleGenesSL = moduleColorsSL==moduleSL;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipSL[moduleGenesSL, column]),
                   abs(geneTraitSignificanceSL[moduleGenesSL, 1]),
                   xlab = paste("Module Membership in", moduleSL, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleSL)

moduleSL = "brown"
column = match(moduleSL, modNamesSL);
moduleGenesSL = moduleColorsSL==moduleSL;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipSL[moduleGenesSL, column]),
                   abs(geneTraitSignificanceSL[moduleGenesSL, 1]),
                   xlab = paste("Module Membership in", moduleSL, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleSL)

moduleSL = "yellow"
column = match(moduleSL, modNamesSL);
moduleGenesSL = moduleColorsSL==moduleSL;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipSL[moduleGenesSL, column]),
                   abs(geneTraitSignificanceSL[moduleGenesSL, 1]),
                   xlab = paste("Module Membership in", moduleSL, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleSL)



library(AnnotationDbi)
library(pathview)
library(GO.db)
library(GOstats)
library(Mus.musculus)
library(org.Mm.eg.db)
library(annotate)

GenesPinkLL <- data.frame(rownames(geneModuleMembershipLL)[moduleColorsLL==moduleLL] %>%
                          mapIds(org.Mm.eg.db, .,  'SYMBOL','ENTREZID',multiVals="first"))

rownames(geneModuleMembershipLL)[moduleColorsLL==moduleLL]


GenesPinkML <- data.frame(rownames(geneModuleMembershipML)[moduleColorsML==moduleML] %>%
                            mapIds(org.Mm.eg.db, .,  'SYMBOL','ENTREZID',multiVals="first"))

rownames(geneModuleMembershipML)[moduleColorsML==moduleML]


GenesPinkSL <- data.frame(rownames(geneModuleMembershipSL)[moduleColorsSL==moduleSL] %>%
                            mapIds(org.Mm.eg.db, .,  'SYMBOL','ENTREZID',multiVals="first"))

rownames(geneModuleMembershipSL)[moduleColorsSL==moduleSL]


# Create the starting data frame
geneInfo0LL = data.frame(geneSymbolLL = rownames(geneModuleMembershipLL)%>%mapIds(org.Mm.eg.db, .,  'SYMBOL','ENTREZID',multiVals="first") ,
                       EntrezIDLL = rownames(geneModuleMembershipLL),
                       moduleColorLL = moduleColorsLL,
                       geneTraitSignificanceLL,
                       GSPvalueLL)


geneInfo0ML = data.frame(geneSymboML = rownames(geneModuleMembershipML)%>%mapIds(org.Mm.eg.db, .,  'SYMBOL','ENTREZID',multiVals="first") ,
                         EntrezIDML = rownames(geneModuleMembershipML),
                         moduleColorML = moduleColorsML,
                         geneTraitSignificanceML,
                         GSPvalueML)


geneInfo0SL = data.frame(geneSymbolSL = rownames(geneModuleMembershipSL)%>%mapIds(org.Mm.eg.db, .,  'SYMBOL','ENTREZID',multiVals="first") ,
                         EntrezIDSL = rownames(geneModuleMembershipSL),
                         moduleColorSL = moduleColorsSL,
                         geneTraitSignificanceSL,
                         GSPvalueSL)



# Order modules by their significance for WAT
modOrderLL = order(-abs(cor(MEsLL, percentageWATLL, use = "p")));

modOrderML = order(-abs(cor(MEsML, percentageWATML, use = "p")));

modOrderSL = order(-abs(cor(MEsSL, percentageWATSL, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembershipLL))
{
  oldNamesLL = names(geneInfo0LL)
  geneInfo0LL = data.frame(geneInfo0LL, geneModuleMembershipLL[, modOrderLL[mod]], 
                         MMPvalueLL[, modOrderLL[mod]]);
  names(geneInfo0LL) = c(oldNamesLL, paste("MM.", modNamesLL[modOrderLL[mod]], sep=""),
                       paste("p.MM.", modNamesLL[modOrderLL[mod]], sep=""))
}


for (mod in 1:ncol(geneModuleMembershipML))
{
  oldNamesML = names(geneInfo0ML)
  geneInfo0ML = data.frame(geneInfo0ML, geneModuleMembershipML[, modOrderML[mod]], 
                           MMPvalueML[, modOrderML[mod]]);
  names(geneInfo0ML) = c(oldNamesML, paste("MM.", modNamesML[modOrderML[mod]], sep=""),
                       paste("p.MM.", modNamesML[modOrderML[mod]], sep=""))
}


for (mod in 1:ncol(geneModuleMembershipSL))
{
  oldNamesSL = names(geneInfo0SL)
  geneInfo0SL = data.frame(geneInfo0SL, geneModuleMembershipSL[, modOrderSL[mod]], 
                           MMPvalueSL[, modOrderSL[mod]]);
  names(geneInfo0SL) = c(oldNamesSL, paste("MM.", modNamesSL[modOrderSL[mod]], sep=""),
                       paste("p.MM.", modNamesSL[modOrderSL[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrderLL = order(geneInfo0LL$moduleColorLL, -abs(geneInfo0LL$GS...WAT.LL));
geneInfo0LL = geneInfo0LL[geneOrderLL, ]

geneOrderML = order(geneInfo0ML$moduleColorML, -abs(geneInfo0ML$GS...WAT.ML));
geneInfo0ML = geneInfo0ML[geneOrderML, ]

geneOrderSL = order(geneInfo0SL$moduleColorSL, -abs(geneInfo0SL$GS...WAT.SL));
geneInfo0SL = geneInfo0SL[geneOrderSL, ]


#
# ===== Code chunk 6 - Save =====
#

write.csv(geneInfo0LL, file = "geneInfo_allGenerationLL.csv")
write.csv(geneInfo0ML, file = "geneInfo_allGenerationML.csv")
write.csv(geneInfo0SL, file = "geneInfo_allGenerationSL.csv")

#write.xlsx(geneInfo, file = "geneInfo_allGeneration.xlsx",sheetName = "AllGen", append = FALSE)

save(geneInfo0LL, file = "03-geneInfo_allGenerationLL.Rdata" )
save(geneInfo0ML, file = "03-geneInfo_allGenerationML.Rdata" )
save(geneInfo0SL, file = "03-geneInfo_allGenerationSL.Rdata" )

# NOTE: 
# GS means gene significance as (the absolute value of) the correlation between the gene and the trait (e.g. WAT).
# MM means module membership. That is a quantitative (fuzzy-based) measure of  membership defined as as the correlation of 
# the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.

