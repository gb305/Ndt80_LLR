####################################################################################################################
#######################  HOTSPOT TABLE GENERATION  ############################################
####################################################################################################################

detach("package:data.table", unload=TRUE) ; detach("package:dplyr", unload=TRUE)
require("e1071") # This pacakge permits smoothing functions (used later)
require(stringr); require(doParallel) ;require(plyr)
options(scipen=999) #Suppresses scientific notation appearing in plots/graphs etc

#Load hotpot template
#Final update template for sae2D ndt80D and sae2D ndt80D tel1D is called 0.193V2_template_UPDATE_3473HS
setwd("ENTER WORKING DIRECTORY")
files = list.files(pattern="template");files
template= read.table(files[2], sep = "\t", header=TRUE) # read template .txt file
template$Length = template$End - template$Start # Define length of hotspots
template$Midpoint = template$Start + (template$Length/2) # Define hotspots'midpoint

# Now point at the data to be processed (.rds format):
setwd("ENTER WORKING DIRECTORY")
parent.directory <- getwd()
files2 <- list.files(pattern = ".rds"); files2

#LOAD DSBL LIST DATA from rds
filesToload = 2L
load(files2[filesToload])
data = DSBList
BG = read.table("BGreads8D.txt", sep = "\t", header=TRUE) #Import background datatable
dflistNames <- as.character(BG$Strain)
BGmean <- as.vector(BG[BG$Strain == dflistNames,"MeanCore"])
nfiles = length(BGmean)

###############################################################################################################################################
############ START HERE ONCE DATAFRAMES ARE LOADED ######################
##########################################################################################################################################################

strains=c(1:nfiles) # Process these numbered dataframes from the DSBList
hotspots=template
extend=300 # To check hotspot and the nearby region (bp)
hs=list() # list of tables containing info on the hotspots for each gene for each strain
cl <- makeCluster(8)
registerDoParallel(cl)
writeLines(c(""), "log.txt")
hs=NULL
r=NULL

r=foreach (k = strains) %dopar%
  { #step through sequentially each dataframe/strain
  DSBList=NULL; Mreads=NULL
  # DSB map files now loaded in within each loop/instance
  DSBList = data[[k]] #Import datatable
  Mreads=sum(DSBList$Watson+DSBList$Crick)/1000000 # Calculate Million reads per sample for conveting to HpM
   DSBList$Watson = DSBList$Watson/Mreads
   DSBList$Crick = DSBList$Crick/Mreads
   
  hs[[k]]=template[c("Chr","Start","End","Length","Midpoint")]
  hs[[k]]=hs[[k]][NULL,]
  row=1
  sink("log.txt", append=TRUE) # Send console output to text file to monitor run

  for (j in 1:16){
    hs1=NULL
    hs1=data.frame(subset(template[1:5], Chr==j)) # subsetting the datatables first by chromosome massively speeds the script up!
    hs1[,c("WatsonHpM","CrickHpM","TotalHpM","BGHpM","Total-BGHpM","WatsonHpM300","CrickHpM300","TotalHpM300","BGHpM300","Total-BGHpM300","NormHpM","NormHpM300","NormHpChr")]=NA #Fill in missing columns before rbind call
    temp1=subset(DSBList, Chr==j) # subsetting the datatables first by chromosome massively speeds the script up! I am sure this is most important for the large DSBList tables!
    hs[[k]]=rbind(hs[[k]],hs1)
    cat("\r", "Job", k, dflistNames[k], "Chromosome", j, "Hotspot", row, "of", nrow(hotspots)); flush.console() # Keep track of progress. cat "\r" overprints to same line of console
    for (i in 1:nrow(hs1)) {
      temp=subset(temp1, Pos>=(hs1[i,"Start"]) & Pos<=(hs1[i,"End"])) # Create temp vector with DSB hits across each hotspot in the table
      hs[[k]][row,"WatsonHpM"]=sum(temp$Watson) # Calculate sum of hits within this region/Mreads[k]
      hs[[k]][row,"CrickHpM"]=sum(temp$Crick) # Calculate sum of hits within this region/Mreads[k]
      
      temp=subset(temp1, Pos>=((hs1[i,"Start"])-extend) & Pos<=((hs1[i,"End"])+extend)) # Create temp vector with DSB hits across hotspot in table +/-300bp
      hs[[k]][row,"WatsonHpM300"]=sum(temp$Watson) # Calculate sum of hits within this region/Mreads[k]
      hs[[k]][row,"CrickHpM300"]=sum(temp$Crick)# Calculate sum of hits within this region/Mreads[k]
  
      row=row+1 # Increment counter
    }
  }
  
  # Vectorised math:
  hs[[k]][1:nrow(hs[[k]]),"TotalHpM"]=hs[[k]][1:nrow(hs[[k]]),"WatsonHpM"]+hs[[k]][1:nrow(hs[[k]]),"CrickHpM"]
  hs[[k]][1:nrow(hs[[k]]),"TotalHpM300"]=hs[[k]][1:nrow(hs[[k]]),"WatsonHpM300"]+hs[[k]][1:nrow(hs[[k]]),"CrickHpM300"]
  hs[[k]][1:nrow(hs[[k]]),"BGHpM"]=BGmean[k]*hs[[k]][1:nrow(hs[[k]]),"Length"]
  hs[[k]][1:nrow(hs[[k]]),"BGHpM300"]=BGmean[k]*(hs[[k]][1:nrow(hs[[k]]),"Length"]+(extend*2))
  hs[[k]][1:nrow(hs[[k]]),"Total-BGHpM"]=hs[[k]][1:nrow(hs[[k]]),"TotalHpM"]-hs[[k]][1:nrow(hs[[k]]),"BGHpM"]
  hs[[k]][1:nrow(hs[[k]]),"Total-BGHpM"][hs[[k]][1:nrow(hs[[k]]),"Total-BGHpM"] < 0] = 0  #NEW - MAKE ALL NEGATIVE VALUES INTO 0
  hs[[k]][1:nrow(hs[[k]]),"Total-BGHpM300"]=hs[[k]][1:nrow(hs[[k]]),"TotalHpM300"]-hs[[k]][1:nrow(hs[[k]]),"BGHpM300"]
  hs[[k]][1:nrow(hs[[k]]),"NormHpM"]=hs[[k]][1:nrow(hs[[k]]),"Total-BGHpM"]
  hs[[k]][1:nrow(hs[[k]]),"NormHpM"][hs[[k]][1:nrow(hs[[k]]),"NormHpM"] < 0] = 0 # Convert all -ve values to zero
  FinalSum=sum(hs[[k]][1:nrow(hs[[k]]),"NormHpM"])
  hs[[k]][1:nrow(hs[[k]]),"NormHpM"]=hs[[k]][1:nrow(hs[[k]]),"NormHpM"]/FinalSum*1000000
  hs[[k]][1:nrow(hs[[k]]),"NormHpM300"]=hs[[k]][1:nrow(hs[[k]]),"Total-BGHpM300"]
  hs[[k]][1:nrow(hs[[k]]),"NormHpM300"][hs[[k]][1:nrow(hs[[k]]),"NormHpM300"] < 0] = 0 # Convert all -ve values to zero
  FinalSum=sum(hs[[k]][1:nrow(hs[[k]]),"NormHpM300"])
  hs[[k]][1:nrow(hs[[k]]),"NormHpM300"]=hs[[k]][1:nrow(hs[[k]]),"NormHpM300"]/FinalSum*1000000
  
  #Loop to calculate NormHpChr
  rowA=0
  for (i in 1:16){
    aChr=subset(hs[[k]],Chr==i)
    aChrDSBs=nrow(aChr) # Number of hotspots per chromosome (rows)
    aChrTotal=sum(aChr[1:aChrDSBs,"NormHpM"])
    hs[[k]][(rowA+1):(aChrDSBs+rowA),"NormHpChr"]=aChr[1:aChrDSBs,"NormHpM"]/aChrTotal*1000000
    rowA=rowA+aChrDSBs
  }
  cat("\r", "Job", k, "COMPLETED", dflistNames[k], "Chromosome", j, "Hotspot", row-1, "of", nrow(hotspots)); flush.console()
  hs[[k]][1:18]=hs[[k]][1:18] # For unknown reasons this code is ESSENTIAL to get the script to populate hs[[k]] with anything. Otherwise it returns "NULL" 
  #hs[[k]][11:23]=round(hs[[k]][11:23], digits=2)
  
  ######### Write tables to text file ##########
  parent.directory
  out = paste("TEST0.193_Hotspot.Table.",dflistNames[k],".txt", sep="") # Set hotspot table name
  write.table(hs[[k]], out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=F)

  }
stopCluster(cl)

sink()
