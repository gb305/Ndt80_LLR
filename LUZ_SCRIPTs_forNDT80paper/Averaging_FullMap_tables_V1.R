#############################################################################
#################### Averaging FullMap tables ###############################
#############################################################################
#Packages to load
require("e1071") ;require(stringr) ; require(doParallel) ; library(doParallel) ; require(plyr)
options(scipen=999) #Suppresses scientific notation appearing in plots/graphs etc

#Set Chromosome size
ChrSize = c(230218,813184,316620+1173+2941 ,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066) # Includes H4L2 insert
ChrSize = ChrSize+1000 #add some padding in case sizes are off
##############################################################################################################################################
# Set working directory to Fulmap files
setwd("ENTER WORKING DIRECTORY")

###############################################################################################################################################
# Import histogram FullMap files for each strain in working directory and tally up the total number of Million mapped reads
Mreads=NULL; dflist=list();dflistNames=NULL
#Read in all tables with string "Full.Map."
files = list.files(pattern="FullMap") ; files # import files names with "FULLMAP" string into variable "files"
strain <- "sae2ndt80tel1" #Change file name
dflistNames = substr(files, 29, nchar(files)-4) # Shorten filename ; ADJUST
dflistNames
nfiles = length(files) # Count number of files
cl <- makeCluster(4)
registerDoParallel(cl)
dflist=foreach (k = 1:nfiles) %dopar% { dflist[[k]] = read.table(files[k], sep = "\t", header=TRUE) } #Import datatable
stopCluster(cl)

###############################################################################################################################################
# Calculate Mreads for each datatable
for (i in 1:nfiles){Mreads[i]=sum(dflist[[i]]$Watson+dflist[[i]]$Crick)/1000000}

# Divide Watson and Crick values by Mreads to convert to HpM
for (i in 1:nfiles){
  dflist[[i]]$Watson=dflist[[i]]$Watson/Mreads[i] # Hash out these lines if you don't want the output converting to HpM 
  dflist[[i]]$Crick=dflist[[i]]$Crick/Mreads[i] # Hash out these lines if you don't want the output converting to HpM 
  }

restored=NULL
output=list()

# Create zeroed out array for entire genome so that FullMaps can be added to it.
for (i in 1:nfiles){
  outputA=NULL
for (chrom in 1:16) { #Step thorugh each chromosome
a2=subset(dflist[[i]], Chr==chrom)
restored <- data.frame(Chr=chrom, Pos=(0:ChrSize[chrom])) # change range if required
restored <- merge(restored,a2, all=TRUE)
restored[is.na(restored)] <- 0
#### NOTE:BE CAREFUL HERE. For RAW data, the following line needs running. But if combining tabels tht are already shifted, it needs hashing out!!!
#restored$Crick=c(tail(restored$Crick,-1),0) ##### Offset the Crick data by -1 bp and fill empty place with zero on the right end of chromosome ### HASH OUT IF COMBINING DUPLICATES ####
#Note: ONLY do preceding line if you want the W and C hits to be on same row - for creating Total column for potential filtering purposes
#Note2: It may be important for the W/C ratio/bias though...
outputA=rbind(outputA,restored) #Bind the chromosomes back together
}
output[[i]]=outputA #Store output in list
}

# Merge files to create single averaged FullMap table
Final=output[[1]]
for (i in 2:nfiles) {
  Final$Watson=Final$Watson+output[[i]]$Watson
  Final$Crick=Final$Crick+output[[i]]$Crick
  }
FinalMreads=sum(Final$Watson+Final$Crick)/1000000  # Check that it adds up correctly

# Convert to HpM
Final$Watson=Final$Watson/FinalMreads # Hash out these lines if you don't want the output converting to HpM 
Final$Crick=Final$Crick/FinalMreads # Hash out these lines if you don't want the output converting to HpM 

#Round values to 4 decimal places
Final[3:4]=round(Final[3:4],4)

#Create output tables with varying thresholds. If a value is below these thresholds for either Watson or Crick it is converted to zero.
#Thresholds=c(0.0, 0.1, 0.2 ,0.3, 0.5, 1.0)
Thresholds=c(0)
FinalTable=list()
for (i in 1:length(Thresholds)) {
  FinalTable[[i]]=subset(Final, Watson >Thresholds[i] | Crick >Thresholds[i])}

#Write out master files EDIT NAME OF STRING!
wd = getwd()
for (i in 1:length(Thresholds)) {
out = paste(strain,"_","CombinedFullMap","_","Average",Thresholds[i],"_",".txt", sep="", collapse="")
write.table(FinalTable[[i]], out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=F)
}

