##########################################################77 SNP#####################################################
setwd("D:/7_Papers/30_T2D/MR/BMI/data_processing/");

library("grid")
library("meta")
library("stringr")

rawData <- read.table('mediaData_leave/78SNPs_005.txt', sep='\t', head=T)          #78SNPs_001.txt   78SNPs_000.txt
#广义逆方差meta分析，基于估计及其标准误差的固定效应和随机效应meta分析
metal <- metagen(rawData[,15], rawData[,16], studlab=rawData[,1], sm="OR");
metal
attributes(metal)
metal$studlab
#漏斗图不对称性检验（对发表偏倚进行定量的检测，p<0.05，证实有发表偏倚）
#k.min进行漏斗图不对称性测试的最少研究次数
#rank秩相关，linreg线性回归
#不知道为啥method写错也能运行
metabias(metal, k.min =1, ethod.bias="linreg")
metabias(metal, k.min =5, ethod.bias="rank")


pdf("D:/1.pdf", width=12, height=12)
funnel(metal)
lines(x=c(1.24,1.24), y=c(1.48, 0), lty="dashed", col="red")
dev.off()

#森林图
#family="GB1",
pdf("D:/3.pdf", width=12, height=12)
forest(metal)
dev.off()

#metal$pval.random
#metal$TE.random
#metal$seTE.random
#metal$data[,2]
#forest$data[,2]

# 
# fileName <- paste('mediaData_leave/', rawData[1,1], sep='');
# fileName <- str_trim(fileName, side='right')
# fileName <- paste(fileName, '.txt', sep='');
# fileName <- str_trim(fileName, side='right')
# str_trim(fileName,side='right')
# 


rawData <- read.table('mediaData_leave/78SNPs_005.txt', sep='\t', head=T)          #78SNPs_001.txt   78SNPs_000.txt
SNPs <- matrix(1:nrow(rawData), ncol=1);
I2 <- matrix(1:nrow(rawData), ncol=1);
P <- matrix(1:nrow(rawData), ncol=1);
OR <- matrix(1:nrow(rawData), ncol=1);
LOR <- matrix(1:nrow(rawData), ncol=1);
UOR <- matrix(1:nrow(rawData), ncol=1);
POR <- matrix(1:nrow(rawData), ncol=1);

#All <- matrix(1:nrow(rawData), ncol=1);

for(i in 1:nrow(rawData)) #
{
  fileName <- paste('mediaData_leave/78_005/', rawData[i,1], sep='');          #78SNPs_001.txt   78SNPs_000.txt
  fileName <- str_trim(fileName, side='right')
  fileName <- paste(fileName, '.txt', sep='');
  fileName <- str_trim(fileName, side='right')
  fileName
  
  oneOut <- read.table(fileName, sep='\t', head=T);
  
  metalOneOut <- metagen(oneOut[,15], oneOut[,16], sm="OR");
  #metalOneOut$pval.fixed;
  #All[i,1] <- metalOneOut
  
  write.table(metacum(metalOneOut), file="mediaData_leave/all.txt", sep="\t",                          #78SNPs_001.txt   78SNPs_000.txt
              append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE);
  #metalOneOut
  #metacum(metalOneOut, pooled="fixed")
  
  SNPs[i, 1] <- str_trim(rawData[i,1]);
  
  P[i, 1] <- 1 - pchisq(metalOneOut$Q, metalOneOut$df.Q, lower.tail = TRUE, log.p = FALSE)
  OR[i, 1] <- exp(metalOneOut$TE.fixed);
  LOR[i, 1] <- exp(metalOneOut$lower.fixed);
  UOR[i, 1] <- exp(metalOneOut$upper.fixed);
  
  
  POR[i, 1] <- metalOneOut$pval.fixed;
  
}

SNPs <- data.frame(SNPs, I2, P, OR, LOR, UOR, POR);

write(c("SNP\tI2\tP\tOR\tLOR\tUOR\tpvalue"), file="mediaData_leave/78_005_oneOut.txt", sep="\t");          #78SNPs_001.txt   78SNPs_000.txt
write.table(SNPs, file="mediaData_leave/78_005_oneOut.txt", sep="\t",                          #78SNPs_001.txt   78SNPs_000.txt
            append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE);

##########################################################77 SNP#####################################################





##########################################################96 SNP#####################################################
setwd("D:/7_Papers/others/??/T2D/data_processing/");

library("grid")
library("meta")
library("stringr")

rawData <- read.table('mediaData/97SNPs.txt', sep='\t', head=T)

metal <- metagen(rawData[,15], rawData[,16], sm="OR");
metal
metal$pval.random
metal$TE.random
metal$seTE.random


fileName <- paste('mediaData/', rawData[1,1], sep='');
fileName <- str_trim(fileName, side='right')
fileName <- paste(fileName, '.txt', sep='');
fileName <- str_trim(fileName, side='right')
str_trim(fileName,side='right')



SNPs <- matrix(1:nrow(rawData), ncol=1);
beta <- matrix(1:nrow(rawData), ncol=1);
se <- matrix(1:nrow(rawData), ncol=1);
randomP <- matrix(1:nrow(rawData), ncol=1);

for(i in 1:nrow(rawData))
{
  fileName <- paste('mediaData/97/', rawData[i,1], sep='');
  fileName <- str_trim(fileName, side='right')
  fileName <- paste(fileName, '.txt', sep='');
  fileName <- str_trim(fileName, side='right')
  fileName
  
  oneOut <- read.table(fileName, sep='\t', head=T);
  
  metalOneOut <- metagen(oneOut[,15], oneOut[,16], sm="OR");
  metalOneOut$pval.random;
  
  SNPs[i, 1] <- str_trim(rawData[i,1]);
  beta[i, 1] <- metalOneOut$TE.random;
  se[i, 1] <- metalOneOut$seTE.random;
  randomP[i, 1] <- metalOneOut$pval.random;
  
}

SNPs <- data.frame(SNPs, beta, se, randomP);

write.table(SNPs, file="mediaData/97_oneOut.txt", sep="\t", 
            append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE);

##########################################################96 SNP#####################################################











##########################################################MendelianRandomization#####################################################
library("MendelianRandomization")
#setwd("D:/7_Papers/others/??/T2D/data_processing/");
#setwd("D:/7_Papers/others/??/17-6-MR/AD/");
setwd("D:/7_Papers/30_T2D/MR/data_processing/");

rawData <- read.table('mediaData_leave/78SNPs_000.txt', sep='\t', head=T)
object <- mr_input(bx=rawData[,7], bxse=rawData[,8], by=rawData[,13], byse=rawData[,14])
mr_allmethods(object)

rawData <- read.table('mediaData_leave/78SNPs_001.txt', sep='\t', head=T)
object <- mr_input(bx=rawData[,7], bxse=rawData[,8], by=rawData[,13], byse=rawData[,14])
mr_allmethods(object)

rawData <- read.table('mediaData_leave/78SNPs_005.txt', sep='\t', head=T)
object <- mr_input(bx=rawData[,7], bxse=rawData[,8], by=rawData[,13], byse=rawData[,14])
mr_allmethods(object)

#pdf("D:/1.pdf", width=12, height=12)
mr_plot(object)
#dev.off()




rawData <- read.table('mediaData/97SNPs_000.txt', sep='\t', head=T)
object <- mr_input(bx=rawData[,7], bxse=rawData[,8], by=rawData[,13], byse=rawData[,14])
mr_allmethods(object)

rawData <- read.table('mediaData/97SNPs_001.txt', sep='\t', head=T)
object <- mr_input(bx=rawData[,7], bxse=rawData[,8], by=rawData[,13], byse=rawData[,14])
mr_allmethods(object)

rawData <- read.table('mediaData/97SNPs_005.txt', sep='\t', head=T)
object <- mr_input(bx=rawData[,7], bxse=rawData[,8], by=rawData[,13], byse=rawData[,14])
mr_allmethods(object)



##########################################################MendelianRandomization#####################################################






























library(grid)
library(meta)
data <- read.table('F://data.txt')
or <- data[,1]
se <- data[,2]
logor <- log(or)
(or.fem <- metagen(logor, se, sm = "OR"))

metabias(or.fem, k.min =1, ethod.bias="linreg")
metabias(or.fem, k.min =5, ethod.bias="rank")
