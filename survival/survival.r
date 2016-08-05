# Copyright (c) 2014-2016  Peter Schüffler, Stefan Bauer [bauers@inf.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
# 
# Survival Analysis with R given the cell classification

########################################################
#Preprocessing
########################################################
require(survival)

X <- read.table("Documents/cellclassification/patches/RCC/ZTMA190505_clean.csv", sep=";", header=TRUE)
sv <- Surv(X$Monate, X$Event)

DI <- read.csv("Desktop/list_files.csv", sep="/",header=TRUE)

library(splitstackshape)

df5<- cSplit(indt = DI, splitCols= "patch.", sep = "_", drop = FALSE)

Net18<- read.csv("Documents/cellclassification/final_classification_18/final_classification_18.csv", header=TRUE)
Net34<- read.csv("Documents/cellclassification/final_classification_34.csv",header=TRUE)

help1<-merge(df5, Net18, by.x = "ID", by.y = "ID18", all = TRUE)
help2<-merge(help1, Net34, by.x = "ID", by.y = "ID34", all = TRUE)

help2$predcancandstain34 <-ifelse((help2$Label34<2) & (help2$patch._5>1) , 34, 0)
help2$predcancandstain18 <-ifelse((help2$Label18<2) & (help2$patch._5>1) , 18, 0)

help2$Label18<-ifelse(help2$Label18>1,0,18) 
help2$Label34<-ifelse(help2$Label34>1,0,34) 

#help2$predcancandstain18 <- as.numeric(as.character(help2$predcancandstain18))

combined18<-aggregate(cbind(help2$predcancandstain18, help2$Label18), list(help2$picture), sum)
combined18$V1<-combined18$V1/18
combined18$V2<-combined18$V2/18

colnames(combined18) <- c("picid18", "stainandcancer18","cancer18")
combined18$propTMA18<-ifelse(combined18$cancer18!=0,combined18$stainandcancer18/combined18$cancer18,0)

combined34<-aggregate(cbind(help2$predcancandstain34,help2$Label34), list(help2$picture), sum)
combined34$V1<-combined34$V1/34
combined34$V2<-combined34$V2/34

colnames(combined34) <-  c("picid34", "stainandcancer34","cancer34")
combined34$propTMA34<-ifelse(combined34$cancer34!=0,combined34$stainandcancer34/combined34$cancer34,0)

help3<-cSplit(indt = X, splitCols  = "Filename", sep = "\\", drop = FALSE)
help4<-cSplit(indt = help3, splitCols  = "Filename_11", sep = ".", drop = TRUE)

help5<-merge(help4, combined18, by.x = "Filename_11_1", by.y = "picid18", all = TRUE)
Y<-merge(help5, combined34, by.x = "Filename_11_1", by.y = "picid34", all = TRUE)

########################################################
#Survival Analysis for Pathologist
###############################################################

########
#Regression: Staining estimation explains the chances for survival
########

coxFit <- coxph(sv~TMA_Mib1_perc, data=X)
waldP <- summary(coxFit)$waldtest[3]

########
#Split (low, (median), high) 
########

###
### two equally sized groups (bin-ary)
###
med <- median(X$TMA_Mib1_perc)
#med <- 2

# P-Value by coxph Wald-test
coxFit_bin <- coxph(sv~ X$TMA_Mib1_perc>=med , data=X)
waldP_bin <- summary(coxFit_bin)$waldtest[3]

# P-Value by survival curve differences (chisquared distributed)

sdf_bin <- survdiff(sv~X$TMA_Mib1_perc>=med) 
p_bin <- 1 - pchisq(sdf_bin$chisq, length(sdf_bin$n) - 1)

########################################################
#Plots 
###############################################################
plot(survfit(sv~X$TMA_Mib1_perc>=med), col=c(1,2), lty=c(1,2), lwd=2,main="Pathologist", xlab="Month", ylab="Cumulative Survival")

leg.txt <- c("1st Quantile", "2nd Quantile")
legend("bottomleft", leg.txt, lty=c(1,2),lwd=2, col=c(1,2), cex=0.9) 

###
# three equally sized groups, based on tertiles
###
X$TMA_Mib1_perc_ter <- cut(X$TMA_Mib1_perc, quantile(X$TMA_Mib1_perc, prob=seq(0,1,1/3)), include.lowest=TRUE, labels=c("1", "2", "3"), right=FALSE)

# P-Value by coxph Wald-test
coxFit_ter <- coxph(sv~ X$TMA_Mib1_perc_ter , data=X)
waldP_ter <- summary(coxFit_ter)$waldtest[3]

# P-Value by survival curve differences (chisquared distributed)
sdf_ter <- survdiff(sv~X$TMA_Mib1_perc_ter) 
p_ter <- 1 - pchisq(sdf_ter$chisq, length(sdf_ter$n) - 1)

# Plot 
plot(survfit(sv~X$TMA_Mib1_perc_ter), col=c(1,2,3), lty=c(1,2,3), lwd=2, main="Pathologist", xlab="Month", ylab="Cumulative Survival")

leg.txt <- c("1st Tertile", "2nd Tertile",
             "3rd Tertile")
legend("bottomleft", leg.txt, lty=c(1,2,3),lwd=2, col=c(1,2,3), cex=0.9) 

########################################################
#Survival Analysis for ResNet34
###############################################################

########
#Regression: Staining explains survival
########
coxFit <- coxph(sv~Y$propTMA34, data=Y)
waldP <- summary(coxFit)$waldtest[3]


########
#Split (low, (median), high) 
########

###
### two equally sized groups (bin-ary)
###

Y$propTMA34<-as.numeric(as.character( Y$propTMA34 ))

med <- median(Y$propTMA34)

# P-Value by coxph Wald-test
coxFit_bin <- coxph(sv~ Y$propTMA34>=med , data=Y)
waldP_bin <- summary(coxFit_bin)$waldtest[3]

# P-Value by survival curve differences (chisquared distributed)

sdf_bin <- survdiff(sv~Y$propTMA34>=med) 
p_bin <- 1 - pchisq(sdf_bin$chisq, length(sdf_bin$n) - 1)

########################################################
#Plots for ResNET34
###############################################################

plot(survfit(sv~Y$propTMA34>=med), col=c(1,2), lty=c(1,2), lwd=2, main="ResNet34", xlab="Month", ylab="Cumulative Survival")

leg.txt <- c("1st Quantile", "2nd Quantile")
legend("bottomleft", leg.txt, lty=c(1,2),lwd=2, col=c(1,2), cex=0.9) 

###
# three equally sized groups, based on tertiles (not really equally sized)
###
Y$propTMA34 <- cut(Y$propTMA34, quantile(Y$propTMA34, prob=seq(0,1,1/3)), include.lowest=TRUE, labels=c("1", "2", "3"), right=FALSE)

# P-Value by coxph Wald-test
coxFit_ter <- coxph(sv~ Y$stained34 , data=Y)
waldP_ter <- summary(coxFit_ter)$waldtest[3]

# P-Value by survival curve differences (chisquared distributed)
sdf_ter <- survdiff(sv~Y$stained34) 
p_ter <- 1 - pchisq(sdf_ter$chisq, length(sdf_ter$n) - 1)

plot(survfit(sv~Y$propTMA34), col=c(1,2,3), lty=c(1,2,3), lwd=2,main="ResNet34", xlab="Month", ylab="Cumulative Survival")

leg.txt <- c("1st Tertile", "2nd Tertile",
             "3rd Tertile")
legend("bottomleft", leg.txt, lty=c(1,2,3),lwd=2, col=c(1,2,3), cex=0.9) 
########################################################
#Survival Analysis for ResNet18
###############################################################

coxFit <- coxph(sv~propTMA18, data=Y)
waldP <- summary(coxFit)$waldtest[3]


########
# Split  (low, (median), high)
########

###
# two equally sized groups (bin-ary)
###

Y$propTMA18<-as.numeric(as.character( Y$propTMA18))

med <- median(Y$propTMA18)


# P-Value by coxph Wald-test
coxFit_bin <- coxph(sv~ Y$propTMA18>=med , data=Y)
waldP_bin <- summary(coxFit_bin)$waldtest[3]

# P-Value by survival curve differences (chisquared distributed)

sdf_bin <- survdiff(sv~Y$propTMA18>=med) 
p_bin <- 1 - pchisq(sdf_bin$chisq, length(sdf_bin$n) - 1)

########################################################
#Plots for ResNET18
###############################################################

# Plot of the two curves, same as in Thomas' Paper
plot(survfit(sv~Y$propTMA18>=med), col=c(1,2), lty=c(1,2), lwd=2,main="ResNet18", xlab="Month", ylab="Cumulative Survival")

leg.txt <- c("1st Quantile", "2nd Quantile")
legend("bottomleft", leg.txt, lty=c(1,2),lwd=2, col=c(1,2), cex=0.9) 


###
# three equally sized groups, based on tertiles 
###
Y$propTMA18 <- cut(Y$propTMA18, quantile(Y$propTMA18, prob=seq(0,1,1/3)), include.lowest=TRUE, labels=c("1", "2", "3"), right=FALSE)

# P-Value by coxph Wald-test
coxFit_ter <- coxph(sv~ Y$propTMA18 , data=Y)
waldP_ter <- summary(coxFit_ter)$waldtest[3]

# P-Value by survival curve differences (chisquared distributed)
sdf_ter <- survdiff(sv~Y$propTMA18) 
p_ter <- 1 - pchisq(sdf_ter$chisq, length(sdf_ter$n) - 1)

plot(survfit(sv~Y$propTMA18), col=c(1,2,3), lty=c(1,2,3), lwd=2, main="ResNet18", xlab="Month", ylab="Cumulative Survival")

leg.txt <- c("1st Tertile", "2nd Tertile",
             "3rd Tertile")
legend("bottomleft", leg.txt, lty=c(1,2,3),lwd=2, col=c(1,2,3), cex=0.9) 