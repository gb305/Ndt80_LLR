########################################################
################NormHpM normalised maps 09.08.19
########################################################
#Load chromosome number and length
setwd("ENTER WORKING DIRECTORY")
chromlengths <- read.table("yeast_chrom_sizes.txt")
chromlengths$kb <- round(chromlengths[,2]/1000, digits=1) # change bp to kb

#Setwt to Hotspot.table file location
setwd("ENTER WORKING DIRECTORY")
parent.directory <- getwd()
setwd(parent.directory); 

# Centromere positions
cen_pos=c(151523.5, 238265, 114443, 449766,152045.5, 148568.5, 496979,105644.5, 355687,
          436366,440187.5, 150887.5, 268090, 628816.5, 326643, 556015) 

##Load data - set each one in a separate variable
files = list.files(pattern="Hotspot.Table") ; files 
DSBListNames = substr(files, 21, nchar(files)-4) ; DSBListNames
#Load in desired order
c1 <- read.table(files[1], sep = "\t", header=TRUE)
c2 <- read.table(files[2], sep = "\t", header=TRUE)
c3 <- read.table(files[3], sep = "\t", header=TRUE)
c4 <- read.table(files[4], sep = "\t", header=TRUE)

contitions <- c("c1","c2","c3","c4") 
#names = DSBListNames
names = c(DSBListNames[1],DSBListNames[2], DSBListNames[3],DSBListNames[34]) # add name in load order 
names(contitions) <- names ; names2 <- names ; contitions 

##Generate a large list that contains all the samples together
DSBList <- list(c1,c2,c3,c4); names(DSBList) <- names  
data <- list(c1,c2,c3,c4); names(data) <- names 

##Samples positions in the data list
sample1 <- 1 ; sample2 <- 2 ; sample3 <- 3 ; sample4 <- 4

#Set strains and order
strains=c(1,2,3,4) 

#Colour vector 
#col.vec <- c("dodgerblue3","goldenrod3","dodgerblue3","goldenrod3") #NEEMAN strains
col.vec <- c("dodgerblue4","dodgerblue1","goldenrod1","goldenrod4") #sae2D ± Ndt80 ±Tel1 strains

#######################################################
##For all chromosomes and strains
########################################################

## Generate Pdf file
pdf(file = "NormHpM.pdf", width = 30, height = 20) 
for(k in 1:16){  
  par(mfrow=c(length(strains),1))
  for(condition in strains){ 
    Chr_coord <- rep(0, chromlengths[k , 2]) #empty vector
    Spo11 <- data.matrix(subset(DSBList[[condition]]$Midpoint, DSBList[[condition]]$Chr == k))
    Spo11_signal <- data.matrix(subset(DSBList[[condition]]$NormHpM, DSBList[[condition]]$Chr == k))
    Chr_coord[Spo11] <- Spo11_signal
    plot(Chr_coord, type = "l", pch=1, lwd = 6, col = col.vec[condition],ylim= c(0,10000), #xlim= c(0,1500000), #Adjust y axis
         main =  paste0(names[condition], "chromosome ", k),xlab = "Chromosome coordinate", ylab = " NormHpM")  
    points(cen_pos[k],0 , col="black", lwd= 3, pch = 19)
  }
}
dev.off()
