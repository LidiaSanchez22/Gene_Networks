#=======================================================================================
 ##Physiological summary
#=======================================================================================

rm(list = ls())

library(xlsx)
library(openxlsx)
library(pacman)
library(readxl)
library(reshape2)
library(ggforce)
library(ggpubr)
library(car)
library(rstatix)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggExtra)
library(magrittr)
library(dplyr)


#datExpr2 <- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx")

#Es necesario haver corrido el Script 1 antes de iniciar el preliminary Data Study
TabPhen<- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx")
datTraits<- TabPhen[, c(1:3, 14:27, 32:34)]

#Convertimos las variables en numéricas dentro del TabPhen (si es necesario)
TabPhen$`INS.SAC`<-as.numeric(TabPhen$`INS.SAC`)

#Dividir en subgrupos
SummaryLL<- TabPhen[TabPhen$Group=="LL",]
SummaryML<- TabPhen[TabPhen$Group=="ML",]
SummarySL<- TabPhen[TabPhen$Group=="SL",]


## DATA SUMMARY
summary(SummaryLL)
summary(SummaryML)
summary(SummarySL)
summary(TabPhen)
summary(datTraits)


#=======================================================================================
 ##Weight vs Time  Plot 
#=======================================================================================

#Weight study

library(xlsx)
library(openxlsx)
library(pacman)
library(readxl)
library(reshape2)
library(ggforce)
library(ggpubr)
library(car)
library(rstatix)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggExtra)
library(magrittr)
library(dplyr)

# Set the default theme to theme_minimal() [in ggplot2]
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

Preliminary <- as.data.frame(read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx"))

# Dataset division into LL, ML, and SL subgroups
LargeLitte <- Preliminary %>% filter(Group == "LL")
MediumLitte <- Preliminary %>% filter(Group == "ML")
SmallLitte <- Preliminary %>% filter(Group == "SL")

# Define custom color palette and prepare the data
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

my3cols <- c("#C20000", "#2E9FDF", "#A36A00",  "#FC4E07","#007500", "#999999", "#FFA500", "#E69F00")
TIME<-c(1,2,3,4,8,12,16,24,32)

##To compute quantile by col
D<-apply(LargeLitte[,5:13],2,quantile,na.rm = T)
##Select 1st and 3rd quartile (25% and 75%) and transpose to make them columns.
D2<-t(D[c(2,4),])
##Median
MED<-apply(LargeLitte[,5:13],2,median, na.rm= T)
Group <- as.character(rep("LargeLitte",9))
DMLL<-as.data.frame(cbind(D2,MED, TIME, Group)) %>% `colnames<-`(c("Q1","Q3","MED","TIME", "Group"))
rm(D,D2,MED,Group)

##To compute quantile by col
D<-apply(MediumLitte[,5:13],2,quantile,na.rm = T)
##Select 1st and 3rd quartile (25% and 75%) and transpose to make them columns.
D2<-t(D[c(2,4),])
##Median
MED<-apply(MediumLitte[,5:13],2,median, na.rm= T)
Group <- as.character(rep("MediumLitte",9))
DMML<-as.data.frame(cbind(D2,MED, TIME, Group)) %>% `colnames<-`(c("Q1","Q3","MED","TIME", "Group"))
rm(D,D2,MED,Group)

##To compute quantile by col
D<-apply(SmallLitte[,5:13],2,quantile,na.rm = T)
##Select 1st and 3rd quartile (25% and 75%) and transpose to make them columns.
D2<-t(D[c(2,4),])
##Median
MED<-apply(SmallLitte[,5:13],2,median, na.rm= T)
Group <- as.character(rep("SmallLitte",9))
DMSL<-as.data.frame(cbind(D2,MED, TIME, Group)) %>% `colnames<-`(c("Q1","Q3","MED","TIME", "Group"))
rm(D,D2,MED,Group)


dfBW <- rbind(DMLL,DMML,DMSL) %>% mutate(FAM = as.factor(Group), TIME=as.numeric(TIME),
                                                          Q1=as.numeric(Q1), Q3=as.numeric(Q3), MED=as.numeric(MED))


##Plot using ggplot2 
pdf("BW_curve.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "US")          # Paper size


ggplot(dfBW, aes(x=TIME, y=MED,group=Group,colour=Group)) + 
  geom_line(size = 1) +
  geom_point(aes(shape =Group), size = 4)+
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.1) +
  labs(x="Time (days)", y = "Body Weight (g)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.2),text = element_text(family = "serif"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 15),
        plot.title = element_text(size = 15)) +
  #scale_fill_manual(values= my3cols)+
  #scale_color_manual(values = my3cols)+ 
  scale_x_continuous(breaks=c(1,2,3,4,8,12,16,24,32), 
  labels=c("0", "7", "14", "21", "30", "60", "90", "120", "180"))


                       
#=================================================================================
  ##Kruskal Wallis Test
#=================================================================================

library(dplyr)
library(devtools)
library(xlsx)
library(openxlsx)
library(pacman)
library(readxl)
library(reshape2)
library(ggforce)
library(ggpubr)
library(car)
library(rstatix)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggExtra)
library(magrittr)
library(dplyr)
library(ggpmisc)
library(ggpubr)

  KWallisTest <- read_excel("C:/Users/lidia/Desktop/LIDIA/MÀSTER EB/ASSIGNATURES/TREBALL FINAL DE MÁSTER/OBESITAT/Table Data New_1.xlsx", sheet="TabData")
    
  head(KWallisTest)
  KWallisTest$Group <- as.factor(KWallisTest$Group)
  levels(KWallisTest$Group)
  
  # Dataset division into LL, ML, and SL subgroups
  LL_Subgroup <- KWallisTest %>% filter(Group == "LL")
  ML_Subgroup <- KWallisTest %>% filter(Group == "ML")
  SL_Subgroup <- KWallisTest %>% filter(Group == "SL")
  
  # BW.SAC----------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(BW.SAC, na.rm = TRUE),
      sd = sd(BW.SAC, na.rm = TRUE),
      median = median(BW.SAC, na.rm = TRUE),
      IQR = IQR(BW.SAC, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "BW.SAC", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "BW.SAC", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "BW.SAC", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "BW.SAC", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(BW.SAC ~ Group, data = KWallisTest)
  
  
  # GLU.SAC---------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(GLU.SAC, na.rm = TRUE),
      sd = sd(GLU.SAC, na.rm = TRUE),
      median = median(GLU.SAC, na.rm = TRUE),
      IQR = IQR(GLU.SAC, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "GLU.SAC", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "GLU.SAC", xlab = "Subgroups", ylim=c(0,100)) 
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
    ggline(KWallisTest, x = "Group", y = "GLU.SAC", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "GLU.SAC", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(GLU.SAC ~ Group, data = KWallisTest)
  
  
  # INS.SAC---------------------------------------------------------------------
  
  KWallisTest$INS.SAC <-as.numeric(KWallisTest$INS.SAC)
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(INS.SAC, na.rm = TRUE),
      sd = sd(INS.SAC, na.rm = TRUE),
      median = median(INS.SAC, na.rm = TRUE),
      IQR = IQR(INS.SAC, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "INS.SAC", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "INS.SAC", xlab = "Subgroups", ylim=c(-1,1))
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "INS.SAC", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "INS.SAC", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(INS.SAC ~ Group, data = KWallisTest)
  
  # INS0MIN---------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(INS0MIN, na.rm = TRUE),
      sd = sd(INS0MIN, na.rm = TRUE),
      median = median(INS0MIN, na.rm = TRUE),
      IQR = IQR(INS0MIN, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "INS0MIN", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "INS0MIN", xlab = "Subgroups", ylim=c(0,1))
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "INS0MIN", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "INS0MIN", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(INS0MIN ~ Group, data = KWallisTest)
  
  # INS30MIN--------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(INS30MIN, na.rm = TRUE),
      sd = sd(INS30MIN, na.rm = TRUE),
      median = median(INS30MIN, na.rm = TRUE),
      IQR = IQR(INS0MIN, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "INS30MIN", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "INS30MIN", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "INS30MIN", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "INS30MIN", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(INS30MIN ~ Group, data = KWallisTest)
  
  # GLU.0-----------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(GLU.0, na.rm = TRUE),
      sd = sd(GLU.0, na.rm = TRUE),
      median = median(GLU.0, na.rm = TRUE),
      IQR = IQR(GLU.0, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "GLU.0", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "GLU.0", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "GLU.0", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "GLU.0", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(GLU.0 ~ Group, data = KWallisTest)
  
  
  # GLU.30----------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(GLU.30, na.rm = TRUE),
      sd = sd(GLU.30, na.rm = TRUE),
      median = median(GLU.0, na.rm = TRUE),
      IQR = IQR(GLU.30, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "GLU.30", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "GLU.30", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "GLU.30", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "GLU.30", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(GLU.30 ~ Group, data = KWallisTest)
  
  
  # D.INS-----------------------------------------------------------------------
  
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(D.INS, na.rm = TRUE),
      sd = sd(D.INS, na.rm = TRUE),
      median = median(D.INS, na.rm = TRUE),
      IQR = IQR(D.INS, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "D.INS", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "D.INS", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "D.INS", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "D.INS", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(D.INS ~ Group, data = KWallisTest)
  
  
  # D.GLU-----------------------------------------------------------------------
  
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(D.GLU, na.rm = TRUE),
      sd = sd(D.GLU, na.rm = TRUE),
      median = median(D.GLU, na.rm = TRUE),
      IQR = IQR(D.GLU, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "D.GLU", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "D.GLU", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "D.GLU", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "D.GLU", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(D.GLU ~ Group, data = KWallisTest)
  
  
  # LIVER-----------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(LIVER, na.rm = TRUE),
      sd = sd(LIVER, na.rm = TRUE),
      median = median(LIVER, na.rm = TRUE),
      IQR = IQR(LIVER, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "LIVER", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "LIVER", xlab = "Subgroups", ylim=c(0,2.5))
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "LIVER", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "LIVER", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(LIVER ~ Group, data = KWallisTest)
  
  
  # WAT-------------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(WAT, na.rm = TRUE),
      sd = sd(WAT, na.rm = TRUE),
      median = median(WAT, na.rm = TRUE),
      IQR = IQR(WAT, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "WAT", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "WAT", xlab = "Subgroups", ylim=c(0,4))
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "WAT", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "WAT", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(WAT ~ Group, data = KWallisTest)
  
  
  # percentage.liver (%.liver)--------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(percentage.liver, na.rm = TRUE),
      sd = sd(percentage.liver, na.rm = TRUE),
      median = median(percentage.liver, na.rm = TRUE),
      IQR = IQR(percentage.liver, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "percentage.liver", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "percentage.liver", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "percentage.liver", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "percentage.liver", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(percentage.liver ~ Group, data = KWallisTest)
  
  # percentage.WAT (% WAT)------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(percentage.WAT, na.rm = TRUE),
      sd = sd(percentage.WAT, na.rm = TRUE),
      median = median(percentage.WAT, na.rm = TRUE),
      IQR = IQR(percentage.WAT, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "percentage.WAT", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "percentage.WAT", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "percentage.WAT", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "percentage.WAT", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(percentage.WAT ~ Group, data = KWallisTest)
  
  
  # TAG.2M----------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(TAG.2M, na.rm = TRUE),
      sd = sd(TAG.2M, na.rm = TRUE),
      median = median(TAG.2M, na.rm = TRUE),
      IQR = IQR(TAG.2M, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "TAG.2M", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "TAG.2M", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "TAG.2M", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "TAG.2M", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(TAG.2M ~ Group, data = KWallisTest)
  
  
  # TAG.4M----------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(TAG.4M, na.rm = TRUE),
      sd = sd(TAG.4M, na.rm = TRUE),
      median = median(TAG.4M, na.rm = TRUE),
      IQR = IQR(TAG.4M, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "TAG.4M", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "TAG.4M", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "TAG.2M", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "TAG.4M", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(TAG.4M ~ Group, data = KWallisTest)
  
  
  # TAG.SAC---------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(TAG.SAC, na.rm = TRUE),
      sd = sd(TAG.SAC, na.rm = TRUE),
      median = median(TAG.SAC, na.rm = TRUE),
      IQR = IQR(TAG.SAC, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "TAG.SAC", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "TAG.SAC", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "TAG.SAC", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "TAG.SAC", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(TAG.SAC ~ Group, data = KWallisTest)
  
  
  # TAG.LivSac------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(TAG.LivSac, na.rm = TRUE),
      sd = sd(TAG.LivSac, na.rm = TRUE),
      median = median(TAG.LivSac, na.rm = TRUE),
      IQR = IQR(TAG.LivSac, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "TAG.LivSac", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "TAG.LivSac", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "TAG.LivSac", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "TAG.LivSac", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(TAG.LivSac ~ Group, data = KWallisTest)
  
  
  # Leptina---------------------------------------------------------------------
 
    group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(Leptina, na.rm = TRUE),
      sd = sd(Leptina, na.rm = TRUE),
      median = median(Leptina, na.rm = TRUE),
      IQR = IQR(Leptina, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "Leptina", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "Leptina", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
    ggline(KWallisTest, x = "Group", y = "Leptina", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "Leptina", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(Leptina ~ Group, data = KWallisTest)
  
  # HOMA.IR---------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(HOMA.IR, na.rm = TRUE),
      sd = sd(HOMA.IR, na.rm = TRUE),
      median = median(HOMA.IR, na.rm = TRUE),
      IQR = IQR(HOMA.IR, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "HOMA.IR", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "HOMA.IR", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
   ggline(KWallisTest, x = "Group", y = "HOMA.IR", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "HOMA.IR", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(HOMA.IR ~ Group, data = KWallisTest)
  
  
  # AUC-------------------------------------------------------------------------
  
  group_by(KWallisTest, Group) %>%
    summarise(
      count = n(),
      mean = mean(AUC, na.rm = TRUE),
      sd = sd(AUC, na.rm = TRUE),
      median = median(AUC, na.rm = TRUE),
      IQR = IQR(AUC, na.rm = TRUE)
    )
  
  # Box plots
  # Plot weight by group and color by group
  ggboxplot(KWallisTest, x = "Group", y = "AUC", 
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("LL", "ML", "SL"),
            ylab = "AUC", xlab = "Subgroups")
  
  # Mean plots
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)
  ggline(KWallisTest, x = "Group", y = "AUC", 
         add = c("mean_se", "jitter"), 
         order = c("LL", "ML", "SL"),
         ylab = "AUC", xlab = "Group")
  
  #Kruskal Wallis Test
  kruskal.test(AUC ~ Group, data = KWallisTest)