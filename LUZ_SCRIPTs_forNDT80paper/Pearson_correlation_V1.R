#############################################################################
# #############################################################################
#Load
library(data.table) ; library(e1071) ; library(doParallel)

#Read hotspot table
Mreads = NULL ; DSBList = list() ; DSBListNames = NULL
cl <- makeCluster(7) ; registerDoParallel(cl)
#Set working directory to where the Hotspot.table files are

setwd("ENTER WORKING DIRECTORY")
files = list.files(pattern = "Hotspot.Table") ; files# import files names with "Hotapot.Table"
DSBListNames = substr(files, 21, nchar(files) - 4) # Shorten file name by x characters from beginning and x characters form end
DSBListNames #Check
nfiles = length(files) # Count number of files
DSBList = foreach (k = 1:nfiles) %dopar% {
  DSBList[[k]] = read.table(files[k], sep = "\t", header = TRUE)
} 
#Import datatable
stopCluster(cl)

#####Pearson correlation

DSBListNames 
rep1 = 1 # Set variable 1
rep2 = 2 # Set variable 2

###NormHpM
pdf(file = paste0(DSBListNames[rep1],"_vs_",DSBListNames[rep2],"_NormHpMANDNormHpChr.pdf"), width = 6, height = 6)
par(mar=c(5,5,5,5)) # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)
plot(DSBList[[rep1]]$NormHpM, 
     DSBList[[rep2]]$NormHpM, type = "p", ylim=c(0.5,22026), log = "xy", 
     xlim =c(0.5,22026),  pch=20, cex = 0.5, col =rgb(0,0,0,0.2),
     xlab=DSBListNames[[rep1]],
     ylab=DSBListNames[[rep2]])
test = lm(DSBList[[rep1]]$NormHpM ~ DSBList[[rep2]]$NormHpM)
title(main = paste0("NormHpM R = ", round(cor(DSBList[[rep1]]$NormHpM, DSBList[[rep2]]$NormHpM, method = "pearson"), 2)," & r2 = ",round(summary(test)$r.squared,2)))
lines(0:10000, 0:10000)
dev.off()

