#SET a main folder with two folders inside fore each mutant
#Add hotspot.table file (without anything in front of hotspot table) inside each mutant folder
#copy main directory and both mutant folder pathways.

rm(list=ls(all=TRUE)) # removes existing data

library(xfun)
# install_github("gb305/nealeLabData",auth_token="63a7cc92827add9796abb9894062ed48d936c979",force=FALSE,ref=)# https://github.com/settings/tokens
library(nealeLabData)
#install_github("gb305/Rpackages",auth_token="63a7cc92827add9796abb9894062ed48d936c979",force=FALSE,ref=)# https://github.com/settings/tokens
library(nealeLabFunctions)

#OTHER PACKAGES:
library(GenomicRanges); library(stats4); library(BiocGenerics); library(S4Vectors)
library(GenomeInfoDb); library(ggplot2) ;library(ggpubr) ;library(GenomicFeatures)
library(AnnotationDbi); library(Biobase) ; library(doParallel) ;library(foreach)
library(iterators) ; library(parallel)

#Set nested folder
setwd("ENTER WORKING DIRECTORY") 

#TOP mutant
s.dir="ENTER WORKING DIRECTORY" 

#BOTTOM mutant
c.dir="ENTER WORKING DIRECTORY" 


CCfoldHS.calc(
  c.dir=c.dir, #BOTTOM
  s.dir=s.dir, #TOP
  p.dir=getwd(),
  coln="NormHpM",
  constant = 5 #Add a constant value to NormHpM before doing the ratio
)

setwd(s.dir) #A .txt file containging log2 FC would be generated inside "TOP" mutant's folder 

files <- list.files(pattern = "NormHpM")
packages.check("doParallel")
packages.check("data.table")
cl <- makeCluster(8)
registerDoParallel(cl)

foreach (i = 1:length(files),.packages=c("nealeLabFunctions","data.table")) %dopar% {

  CCfoldHS.smoother(files[i])}

stopCluster(cl)
setwd(paste0(getwd(),"/smoothed"))
files <- list.files(pattern = ".Rbin")

packages.check("doParallel")
cl <- makeCluster(7)
registerDoParallel(cl)
packages.check("GenomicRanges",BiocM=T)
packages.check("ggplot2")
packages.check("ggpubr")
packages.check("GenomicFeatures",BiocM=T)
 

 foreach (i = 1:length(files),.packages=c("nealeLabFunctions","GenomicRanges","ggplot2","ggpubr")) %dopar% {
CCfoldHS.binning(files,i,binwidth.v=c(100,50,25,10,5), plotmode = F)
  print(i)
 }
stopCluster(cl)

#UP to here a "Smoothed" folder would be generated, inside BinData and BinPlots folders with 5,10,25,50,100kb folders'

##This line generates a CENpileup folder with the.png file insithe the binwith folder you specify
#colhigh =  , collow = arguments to change colours
CCfoldHS.binpileup(log2FoldMax=2,binwidth.v=c(50),invert.col=F,fixmax=T, colhigh = "dodgerblue3", collow = "goldenrod2") 


# TO PLOT REC114:
 setwd("ENTER WORKING DIRECTORY")
 CCfoldHS.binpileup(log2FoldMax=0.3,binwidth.v=c(50),invert.col=T,fixmax=T, 
                    log2foldoffset = 2.5,colhigh = "firebrick", collow = "deepskyblue4")
 
