# AIM : Interfacing network analysis with other data such as functional annotation and gene ontology

rm(list = ls())
#options(java.parameters = "-Xmx1024m")
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(xlsx)
library(readxl)
library(ggvenn)


#
# ===== Code chunk 1 - Load Data =====
#

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/Documents/JIMENEZJOS_02/Res_WGCNA";
setwd(workingDir); 
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "C57-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames2 = load(file = "Mouse-networkConstruction-auto.RData");
lnames2
# Load Modules data saved in the 3rd part.
lnamesLL = load(file = "03-geneInfo_allGenerationLL.Rdata");
lnamesML = load(file = "03-geneInfo_allGenerationML.Rdata");
lnamesSL = load(file = "03-geneInfo_allGenerationSL.Rdata");
lnamesLL
lnamesML
lnamesSL
rm(lnamesLL, lnamesML, lnamesSL, lnames, lnames2)

#
# ==== Code chunk 2 - Output gene lists for use with online software and services ====
#

# To export a list of gene identifiers that can be used as input for several popular gene ontology
# and functional enrichment analysis suites such as David or AmiGO. For example, we write out the EntrezID
# codes for the module highly correlated to WAT into a file.


# Get the corresponding Locuis Link IDs
allLLIDsLL = geneInfo0LL$EntrezIDLL;
allLLIDsML = geneInfo0ML$EntrezIDML;
allLLIDsSL = geneInfo0SL$EntrezIDSL;

# $ Choose interesting modules
intModulesLL = c("darkmagenta", "pink")
for (module in intModulesLL)
{
  # Select module probes
  modGenesLL = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDsLL = allLLIDsLL[modGenesLL];
  # Write them into a file
  fileName = paste("LocusLinkIDsLL-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDsLL), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

intModulesML = c("pink", "black")
for (module in intModulesML)
{
  # Select module probes
  modGenesML = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDsML = allLLIDsML[modGenesML];
  # Write them into a file
  fileName = paste("LocusLinkIDsML-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDsML), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

intModulesSL = c("brown")
for (module in intModulesSL)
{
  # Select module probes
  modGenesSL = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDsSL = allLLIDsSL[modGenesLL];
  # Write them into a file
  fileName = paste("LocusLinkIDsSL-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDsSL), file = fileName,
              row.names = FALSE, col.names = FALSE)
}


# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-allLL.txt", sep="");
write.table(as.data.frame(allLLIDsLL), file = fileName,
            row.names = FALSE, col.names = FALSE)

fileName = paste("LocusLinkIDs-allML.txt", sep="");
write.table(as.data.frame(allLLIDsML), file = fileName,
            row.names = FALSE, col.names = FALSE)

fileName = paste("LocusLinkIDs-allSL.txt", sep="");
write.table(as.data.frame(allLLIDsSL), file = fileName,
            row.names = FALSE, col.names = FALSE)

#
# ===== Code chunk 3 - Enrichment analysis directly within R =====
#

library(AnnotationDbi)
library(pathview)
library(GO.db)
library(GOstats)
library(Mus.musculus)
library(org.Mm.eg.db)
library(annotate)


GOenrLL = GOenrichmentAnalysis(moduleColors, allLLIDsLL, organism = "mouse", nBestP = 10);
GOenrML = GOenrichmentAnalysis(moduleColors, allLLIDsML, organism = "mouse", nBestP = 10);
GOenrSL = GOenrichmentAnalysis(moduleColors, allLLIDsSL, organism = "mouse", nBestP = 10);

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


tabLL = GOenrLL$bestPTerms[[4]]$enrichment
tabML = GOenrML$bestPTerms[[4]]$enrichment
tabSL = GOenrSL$bestPTerms[[4]]$enrichment

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


names(tabLL)
names(tabML)
names(tabSL)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


#write.table(tabLL, file = "GOEnrichmentTableLL.csv", sep = ",", quote = TRUE, row.names = FALSE)
#write.xlsx(geneInfo0LL, file = "GOEnrichmentTableLL.xlsx",sheetName = "AllGenLL", append = FALSE)

save(geneInfo0LL, file = "GoEnrichmentTableLL.RData")


#write.table(tabML, file = "GOEnrichmentTableML.csv", sep = ",", quote = TRUE, row.names = FALSE)
#write.xlsx(geneInfo0ML, file = "GOEnrichmentTableML.xlsx",sheetName = "AllGenML", append = FALSE)

save(geneInfo0LL, file = "GoEnrichmentTableML.RData")



#write.table(tabSL, file = "GOEnrichmentTableSL.csv", sep = ",", quote = TRUE, row.names = FALSE)
#write.xlsx(geneInfo0SL, file = "GOEnrichmentTableSL.xlsx",sheetName = "AllGenSL", append = FALSE)

save(geneInfo0LL, file = "GoEnrichmentTableSL.RData")

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


keepColsLL = c(1, 2, 5, 6, 7, 12, 13);
screenTabLL = tabLL[, keepColsLL];
# Round the numeric columns to 2 decimal places:
numColsLL = c(3, 4);
screenTabLL[, numColsLL] = signif(apply(screenTabLL[, numColsLL], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTabLL[, 7] = substring(screenTabLL[, 7], 1, 40)
# Shorten the column names:
colnames(screenTabLL) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTabLL) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(max.print=1000000)
options(width=1000)
# Finally, display the enrichment table:
screenTabLL


keepColsML = c(1, 2, 5, 6, 7, 12, 13);
screenTabML = tabML[, keepColsML];
# Round the numeric columns to 2 decimal places:
numColsML = c(3, 4);
screenTabML[, numColsML] = signif(apply(screenTabML[, numColsML], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTabML[, 7] = substring(screenTabML[, 7], 1, 40)
# Shorten the column names:
colnames(screenTabML) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTabML) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(max.print=1000000)
options(width=1000)
# Finally, display the enrichment table:
screenTabML


keepColsSL = c(1, 2, 5, 6, 7, 12, 13);
screenTabSL = tabSL[, keepColsSL];
# Round the numeric columns to 2 decimal places:
numColsSL = c(3, 4);
screenTabSL[, numColsSL] = signif(apply(screenTabSL[, numColsSL], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTabSL[, 7] = substring(screenTabSL[, 7], 1, 40)
# Shorten the column names:
colnames(screenTabSL) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTabSL) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(max.print=1000000)
options(width=1000)
# Finally, display the enrichment table:
screenTabSL


#                                                                          #
# ===== Code chunk 7 - Venn Diagram y Enrichment Analysis GO and KEGG =====#
#                                                                          #

#Charging all libraries

library(org.Mm.eg.db)
library(ggvenn)
library(clusterProfiler)
library(airway)
library(edgeR)
library(enrichplot)
library(VennDiagram)
library(AnnotationDbi)
library(cowplot)

## Sacar las muestras y el número de genes de los 3 dataframes
#Combination 1

subgeneInfo0LLa<- geneInfo0LL[geneInfo0LL$moduleColor=="darkmagenta",]
subgeneInfo0MLa<- geneInfo0ML[geneInfo0ML$moduleColor=="black",]
subgeneInfo0SL<- geneInfo0SL[geneInfo0SL$moduleColor== "brown",]

dim(subgeneInfo0LLa)
dim(subgeneInfo0MLa)
dim(subgeneInfo0SL)

#Combination 2

subgeneInfo0LLb<- geneInfo0LL[geneInfo0LL$moduleColor=="pink",]
subgeneInfo0MLb<- geneInfo0ML[geneInfo0ML$moduleColor=="black",]
subgeneInfo0SL<- geneInfo0SL[geneInfo0SL$moduleColor== "brown",]

dim(subgeneInfo0LLb)
dim(subgeneInfo0MLb)
dim(subgeneInfo0SL)

#Combination 3

subgeneInfo0LLc<- geneInfo0LL[geneInfo0LL$moduleColor=="darkmagenta",]
subgeneInfo0MLc<- geneInfo0ML[geneInfo0ML$moduleColor=="pink",]
subgeneInfo0SL<- geneInfo0SL[geneInfo0SL$moduleColor== "brown",]

dim(subgeneInfo0LLc)
dim(subgeneInfo0MLc)
dim(subgeneInfo0SL)

#Combination 4

subgeneInfo0LLd<- geneInfo0LL[geneInfo0LL$moduleColor=="pink",]
subgeneInfo0MLd<- geneInfo0ML[geneInfo0ML$moduleColor=="pink",]
subgeneInfo0SL<- geneInfo0SL[geneInfo0SL$moduleColor== "brown",]

dim(subgeneInfo0LLd)
dim(subgeneInfo0MLd)
dim(subgeneInfo0SL)


--------------------------------------------------------------------------------
  ## VennDiagram ##
  --------------------------------------------------------------------------------
 #Combination 1  
  genes <- paste("gene",1:1000,sep="")
Venn <- list(
  LL= subgeneInfo0LLa$EntrezIDLL,
  ML= subgeneInfo0MLa$EntrezIDML,
  SL= subgeneInfo0SL$EntrezIDSL)

ggvenn(
  Venn, 
  fill_color = c("darkmagenta", "black", "brown"),
  stroke_size = 0.5, set_name_size = 4
)

#Extract the genes from the intersection
Intersect1 <- function (Venn) {  
  # Multiple set version of intersect
  # x is a list
  if (length(Venn) == 1) {
    unlist(Venn)
  } else if (length(Venn) == 2) {
    intersect(Venn[[1]], Venn[[2]])
  } else if (length(Venn) > 2){
    intersect(Venn[[1]], Intersect(Venn[-1]))
  }
}

xx.1 <- list( LL= subgeneInfo0LLa$EntrezIDLL,
              ML= subgeneInfo0MLa$EntrezIDML,
              SL= subgeneInfo0SL$EntrezIDSL)
Intersect1(xx.1)


#Combination 2  
genes <- paste("gene",1:1000,sep="")
Venn <- list(
  LL= subgeneInfo0LLb$EntrezIDLL,
  ML= subgeneInfo0MLb$EntrezIDML,
  SL= subgeneInfo0SL$EntrezIDSL)

ggvenn(
  Venn, 
  fill_color = c("pink", "black", "brown"),
  stroke_size = 0.5, set_name_size = 4
)

#Extract the genes from the intercection
Intersect2 <- function (Venn) {  
  # Multiple set version of intersect
  # x is a list
  if (length(Venn) == 1) {
    unlist(Venn)
  } else if (length(Venn) == 2) {
    intersect(Venn[[1]], Venn[[2]])
  } else if (length(Venn) > 2){
    intersect(Venn[[1]], intersect(Venn[-1]))
  }
}

xx.2 <- list( LL= subgeneInfo0LLb$EntrezIDLL,
              ML= subgeneInfo0MLb$EntrezIDML,
              SL= subgeneInfo0SL$EntrezIDSL)
Intersect2(xx.2)

#Combination 3  
genes <- paste("gene",1:1000,sep="")
Venn <- list(
  LL= subgeneInfo0LLc$EntrezIDLL,
  ML= subgeneInfo0MLc$EntrezIDML,
  SL= subgeneInfo0SL$EntrezIDSL)

ggvenn(
  Venn, 
  fill_color = c("darkmagenta", "pink", "brown"),
  stroke_size = 0.5, set_name_size = 4
)

#Extract the genes from the intercection
Intersect3 <- function (Venn) {  
  # Multiple set version of intersect
  # x is a list
  if (length(Venn) == 1) {
    unlist(Venn)
  } else if (length(Venn) == 2) {
    intersect(Venn[[1]], Venn[[2]])
  } else if (length(Venn) > 2){
    intersect(Venn[[1]], Intersect(Venn[-1]))
  }
}

xx.3 <- list( LL= subgeneInfo0LLc$EntrezIDLL,
              ML= subgeneInfo0MLc$EntrezIDML,
              SL= subgeneInfo0SL$EntrezIDSL)
Intersect3(xx.3)

#Combination 4  
genes <- paste("gene",1:1000,sep="")
Venn <- list(
  LL= subgeneInfo0LLd$EntrezIDLL,
  ML= subgeneInfo0MLd$EntrezIDML,
  SL= subgeneInfo0SL$EntrezIDSL)

ggvenn(
  Venn, 
  fill_color = c("pink", "pink", "brown"),
  stroke_size = 0.5, set_name_size = 4
)

#Extract the genes from the intercection
Intersect4 <- function (Venn) {  
  # Multiple set version of intersect
  # x is a list
  if (length(Venn) == 1) {
    unlist(Venn)
  } else if (length(Venn) == 2) {
    intersect(Venn[[1]], Venn[[2]])
  } else if (length(Venn) > 2){
    intersect(Venn[[1]], Intersect(Venn[-1]))
  }
}

xx.4 <- list( LL= subgeneInfo0LLd$EntrezIDLL,
              ML= subgeneInfo0MLd$EntrezIDML,
              SL= subgeneInfo0SL$EntrezIDSL)
Intersect4(xx.4)

--------------------------------------------------------------------------------
  ## Enrichment Analysis##
  --------------------------------------------------------------------------------
  
  ## GO term Analysis ##
  
  # Pou2f1 is 18986
  # Camk2n2 is 73047
  # Clcf1 is 56708
  
intersecgenes<- c("Pou2f1", "Camk2n2", "Clcf1")
intersecgenes2<- c("18986", "73047", "56708")
#allgenes<- subgeneInfo0SL$EntrezIDSL
intersecgenesc2<- c("71602")
intersecgenes2<- c("Myo1e")

enrich_result <- enrichGO(
  gene = intersecgenesc2,
  OrgDb = org.Mm.eg.db,  # Base de datos de genes
  keyType = "SYMBOL",    # Tipo de identificador (puede variar según tu base de datos)
  ont = "BP",            # Ontología GO (por ejemplo, BP para procesos biológicos)
  pAdjustMethod = "BH",  # Método de ajuste para valores p
  qvalueCutoff = 0.05    # Umbral para el valor q (valor ajustado)
)

tab.go <- as.data.frame(enrich_result)
#tab.go <- subset(tab.go, Count>5) Si lo aplico me queda tabla Nula
tab.go[1:5, 1:6]

head(enrich_result)

#BP stands for Biological Process, MP stands for Molecular Process and CC for Cellular Component
#(The three categories for an GO enrichment analysis) #MusMusculus is Mm


## KEGG analysis##

#mmu is for mouse MusMusculus

enrich_result2 <- enrichKEGG(gene = intersecgenesc2, organism = 'mmu', qvalueCutoff = 0.05)
head(enrich_result2)


tab.kegg <- as.data.frame(enrich_result2)
tab.kegg<- subset(tab.kegg, Count>5)
tab.kegg[1:5, 1:6]



--------------------------------------------------------------------------------
  ## Data Visualization##
--------------------------------------------------------------------------------
  
p1 <- plot(barplot(enrich_result, showCategory = 20))
p1

p2 <- dotplot(enrich_result, showCategory=20) + ggtitle("GO")
p2
p3 <- dotplot(enrich_result, showCategory=20) + ggtitle("Disease")
p3
plot_grid(p2, p3, nrow=2)
p4 <- upsetplot(enrich_result)
p4
p5 <- emapplot(enrich_result)
p5

cowplot::plot_grid(p1, p3, p5, ncol=2, labels=LETTERS[1:3])
