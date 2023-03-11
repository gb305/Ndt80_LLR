###############################################################################################################################################
############ CALCULATING BACKGROUND READS ######################

  #CHANGE TO WORKING DIRECTORY CONTAINING THE INPUT DATA MAPS
  # Load required dataset files:
  setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES")
  AllElementsDUB = read.table("AllElementsDUB_H4L2_Brar_2016.08.16.txt", sep = "\t", header=TRUE) #Import datatable

  # Now point at the data to be processed:
  setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/MAPPING/REC114/sae2Dndt80D background/AVERAGES")
  
  ###############################################################################################################################################
  require("e1071") ; require(stringr);require(doParallel); require(plyr) #Packages to load
  options(scipen=999) #Suppresses scientific notation appearing in plots/graphs etc

  ###############################################################################################################################################
  # Import histogram FullMap files for each strain in working directory and tally up the total number of Million mapped reads
  Mreads=NULL; DSBList=list();DSBListNames=NULL
  #Read in all tables with string "Full.Map."
  files = list.files(pattern="FullMap") ;  files # look for files 
  DSBListNames = c("tel1D", "wt") #Set strain names
  nfiles = length(files) # Count number of files

###############################################################################################################################################
############ START HERE ONCE DATAFRAMES ARE LOADED ######################
##########################################################################################################################################################
# New loop to plot multiple comparisions
strains=c(nfiles) # Analyse these numbered dataframes from the dflist
strains=1:nfiles
###############################################################################################################################################
# MODULE for pulling out largest genes
AllElementsDUB$length=AllElementsDUB$stop-AllElementsDUB$start # Add length of feature column
genes=AllElementsDUB #First make a copy of the ALLElements table
genes=subset(genes, type=="gene" & genename!="Dubious_ORF" & chr !="chrmt")
genes=subset(genes, !genes$genename %in% c("TEL1","NUM1","YRF1-7","YRF1-6","YRF1-3","URA2","TOR2")) # Exclude these genes (emprically determined to be outliers)
genes=subset(genes, length>=5500)
at=sum(genes$length)
genes$start2=genes$start+1000
genes$stop2=genes$stop-1000
genes$length2=genes$stop2-genes$start2
at2=sum(genes$length2)

###############################################################################################################################################

bg=list() # list of tables containing info on the backgrounds for each gene for each strain
DSBList=list() # list of tables of DSB maps
Mreads1=list()
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

r=NULL
r=foreach (k = strains) %dopar% { #step through sequentially each dataframe/strain. r collects within it a list of the foreach loops
  # DSB map files now loaded in within each loop/instance
  DSBList = read.table(files[k], sep = "\t", header=TRUE) #Import datatable
  Mreads=sum(DSBList$Watson+DSBList$Crick)/1000000 # Calculate Million reads per sample for conveting to HpM
  #Mreads1[[k]][1]=Mreads
  
  bg[[k]]=data.frame(NULL)
  bg[[k]][1:nrow(genes),"Genename"]=genes$genename
  #sink("log.txt", append=TRUE) # Send console output to text file to monitor run
  #print(dflistNames[k])
  
for (i in 1:nrow(genes)) {
  cat("\r", i, "of", nrow(genes)); flush.console() # Keep track of progress. cat "\r" overprints to same line of console
  
  temp=subset(DSBList, Chr==(genes[i,"chr"]) & Pos>=(genes[i,"start"]) & Pos<=(genes[i,"stop"])) # Create temp vector with DSB hits across each gene in the table
  bg[[k]][i,"Total"]=sum(temp$Watson+temp$Crick) # Calculate sum of hits within this region
  bg[[k]][i,"Density"]=bg[[k]][i,"Total"]/Mreads/genes[i,"length"] # Calculate density/bp of hits within this region
  temp2=subset(DSBList, Chr==(genes[i,"chr"]) & Pos>=(genes[i,"start2"]) & Pos<=(genes[i,"stop2"])) # Create temp vector with DSB hits across each gene in the table for "core" region
  bg[[k]][i,"TotalCore"]=sum(temp2$Watson+temp2$Crick) # Calculate sum of hits within this region
  bg[[k]][i,"DensityCore"]=bg[[k]][i,"TotalCore"]/Mreads/genes[i,"length2"] # Calculate density/bp of hits within this region
}
  bg[[k]][1:5]=bg[[k]][1:5] # This line is essential in multicore loops for some reason

}
bg=r # Collect the foreach loop
stopCluster(cl)

######################################################################
########## Drawing DotPlots of background counts ####################
######################################################################

wd = getwd(); out = paste(wd,"/","BackgroundReads",Sys.time(),".pdf",sep=""); pdf(file=out, width=15,height=9);
layout(matrix(c(1,2), 1, 2, byrow = T))
for (k in strains) { #step through sequentially each dataframe/strain
dotchart(bg[[k]]$Density, labels=bg[[k]]$Genename, main=DSBListNames[k], xlim=c(0,max(bg[[k]]$Density)), xlab="Total hits per Million reads per bp")
dotchart(bg[[k]]$DensityCore, labels=bg[[k]]$Genename, main=paste(c(DSBListNames[k],"Core")), xlim=c(0,max(bg[[k]]$DensityCore)), xlab="Total hits per Million reads per bp")
}
dev.off()

######################################################################
########## Making background count table "BGreads" ####################
######################################################################

BG=data.frame(NULL)
for (k in strains) {
BG[k,"Strain"]=DSBListNames[k]
BG[k,"Mean"]=mean(bg[[k]]$Density)
BG[k,"StDev"]=sd(bg[[k]]$Density)
BG[k,"StDev%"]=sd(bg[[k]]$Density)/BG[k,"Mean"]*100
BG[k,"MeanCore"]=mean(bg[[k]]$DensityCore)
BG[k,"StDevCore"]=sd(bg[[k]]$DensityCore)
BG[k,"StDevCore%"]=sd(bg[[k]]$DensityCore)/BG[k,"MeanCore"]*100
}
BG[2:7]=signif(BG[2:7], digits=4)

#Write out master files
wd = getwd()
out = paste(wd,"/","BGreads.txt", sep="")
write.table(BG, out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=F)




