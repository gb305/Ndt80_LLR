#########################################################
##################Plot Total HpM hotspot#################
#########################################################

###############################################################################################################################################
#Packages to load
library("e1071") ; library(stringr) ; library(doParallel) ; library(plyr)
require(stringr); require("e1071") # This pacakge permits smoothing functions (used later)
options(scipen=999) #Suppresses scientific notation appearing in plots/graphs etc
require(doParallel) ; require(plyr) 

###############################################################################################################################################
# Load required dataset files:
setwd("ENTER WORKING DIRECTORY")
AllElementsDUB = read.table("AllElementsDUB_H4L2_Brar_2016.08.16.txt", sep = "\t", header=TRUE) #Import datatable
nuc=read.table("nucleoChemUnique_H4L2_2016.08.11.txt", sep = "\t", header=TRUE) #Import datatable
rmm=read.table("RMMSubSamp_simple_H4L2_2016.08.11.txt", sep = "\t", header=TRUE) #Import datatable
chromlengths <- read.table("yeast_chrom_sizes.txt")
rmm1=rmm[seq(1, NROW(rmm), by = 10),] # Subsample every 10th row to reduce resoltuion of RMM plot

# Centromere positions
CEN=c(151523.5, 238265, 114443, 449766,152045.5, 148568.5, 496979,105644.5, 
      355687, 436366,440187.5, 150887.5, 268090, 628816.5, 326643, 556015) 

# Now point at the data to be processed:
setwd("ENTER WORKING DIRECTORY")
BG = read.table("BGreads8D.txt", sep = "\t", header=TRUE) #Import hotspot datatable

# Creaty empty lists
Mreads=NULL; DSBList=list();DSBListNames=NULL
cl <- makeCluster(4)
registerDoParallel(cl)

#Read in all tables with string "FullMap."
files = list.files(pattern="FullMap") ; files # import files names with a pattern string into variable "Fullmap" 
DSBListNames = as.character(BG$Strain) ; DSBListNames
nfiles = length(files) # Count number of files

#Import datatable
DSBList = foreach (k = 1:nfiles) %dopar% { DSBList[[k]] = read.table(files[k], sep = "\t", header=TRUE) } 
stopCluster(cl)

#Bakcground reads
BGmean=NULL; for (i in 1:nfiles){BGmean[i]=unlist(subset(BG, Strain==DSBListNames[i], MeanCore))} # Ensure that background vector BGmean is using same indexe numbering as hit data 
Nfactor=(1-(BGmean*12.01)) # Normalisation factor: based on number of reads that appear NOT to be background
#Nfactor = rep(1, length(DSBList))
for (i in 1:nfiles){Mreads[i]=sum(DSBList[[i]]$Watson+DSBList[[i]]$Crick)/1000000} # Calculate Million reads per sample for conveting to HpM

###############################################################################################################################################
############ START HERE ONCE DATAFRAMES ARE LOADED ######################
###############################################################################################################################################

# MODULE for pulling out specific locus of interest (gene names inside AllElementsDUB_H4L2_Brar_2016.08.16.txt)
##For His4::LEU'2 locus - "H4L2insert" ; For ARE1 - "ARE1" 

orf="ARE1" #Modify here the locus to plot
genes=AllElementsDUB #First make a copy of the ALLElements table
upstream= 4000; downstream= 4000 # bp to extend by in either direction of ORF ; ADJUST
genes=subset(genes, genename==orf | sysname==orf)
xl1=genes$start-upstream
xl2=genes$stop+downstream
chrom=genes$chr

################
# MODULES for overlaying optional datatracks. Use =1 to plot, and other value not to plot
nucflag=1 ; rmmflag=0 ; bgsubtract=1 #minimum value removed from datatracks during plotting stage (NEW)
##########################################################################################################################################################
# New loop to plot multiple comparisions
DSBListNames

strains=c(1,2) # Plot these numbered dataframes from the DSBList **adjust for order you want to plot in
scalar.u=c(1,1,1,1) # Unique scaling factor for each strain in the DSBList (default =1 is identical scaling) **adjust for number of strains
scalar.u=scalar.u*1/Nfactor # Multiply scaling factor by apparent hit reads

col.vec <- c("dodgerblue1","goldenrod1","green4","firebrick") #Adjust colour vector

###############################################################################################################################################
#Plotting: first set up how the plots are organised. How many panes per image for example using the layout command
plotnumber=length(strains) # Number from 1 to 5
if (plotnumber==1) {layout(matrix(c(1,1,1,2,2),5, 1, byrow = T))}
if (plotnumber==2) {layout(matrix(c(1,1,1,2,2,2,3,3), 8, 1, byrow = T))}
if (plotnumber==3) {layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4),11, 1, byrow = T))}
if (plotnumber==4) {layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),14, 1, byrow = T))}
if (plotnumber==5) {layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6),17, 1, byrow = T))}

par(mar=c(1,5,1,0),oma = c(0, 1, 1, 1),las=1) # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)
layout.show((length(strains)+1))

for (k in strains){ #step through sequentially each dataframe/strain

#Subset for region of interest
sae2.0=subset(DSBList[[k]], Chr==chrom & Pos>=xl1 & Pos <=xl2) #Make a sub-table of the sae2-DSB data that only contains those rows where chr = 1 in range of interest

##########################################################################################################################################################
#Decompression code here
sae2.1 <- data.frame(Chr=chrom, Pos=(xl1:xl2)) # Creates expanded empty dataframe with Chr and Pos locations
sae2.1 <- merge(sae2.1,sae2.0, all=TRUE) # Merge expanded empty dataframe with compressed sae2.1 dataframe
sae2.1[is.na(sae2.1)] <- 0 # Convert all NA values to zero
sae2.1$Total=c(head(sae2.1$Watson,-1)+tail(sae2.1$Crick,-1),0) #Total now is adjusted for the 1bp offset of W and Crick
sae2.1[sae2.1$Watson<=bgsubtract,"Watson"]=0 #Remove all signal 2 or below from plots
sae2.1[sae2.1$Crick<=bgsubtract,"Crick"]=0 #Remove all signal 2 or below from plots
sae2.1[sae2.1$Total<=bgsubtract*2,"Total"]=0 #Remove all signal 2 or below from plots
# NOTE: I think it woudl be better to first correct the W/C offset, then remove all bases form W and C where there are not at least 1 hit on both strands - maybe...

##########################################################################################################################################################
# Smoothing function #### temp and smooth are just two temporary vectors. New version creates two smoothed plots for each profile for overlaying
win=1 # hanning window size [1]
win2=101 # hanning window size [2] for overlay
scalar=15 #scalar [1]
scalar2=2 # scalar [2] for overlay
sae2scalar=scalar.u[k]*scalar # adjust this if needed when adjusting hann window smoothing
sae2scalar2=scalar.u[k]*scalar2 # adjust this if needed when adjusting hann window smoothing
hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
hw2=hanning.window(win2) #create hanning window (require package e1071 to be loaded)

##########################################################################################################################################################
temp=NULL
for (j in 3:5){
temp=c(rep(0,win),sae2.1[1:nrow(sae2.1),j], rep(0,win)) # Create vector length of chromosome and extend by the length of the slidign window with zeros at both ends
smooth=filter(temp,hw) # smooth the temp vector using the hann window
smooth=smooth[(win+1):(length(smooth)-win)] # trim smooth to correct lengthsmooth2=smooth2[(win2+1):(length(smooth2)-win2)] # trim smooth to correct length
sae2.1[j+3]=smooth
}
temp2=NULL
for (j in 3:5){
  temp2=c(rep(0,win2),sae2.1[1:nrow(sae2.1),j], rep(0,win2)) # Create vector length of chromosome and extend by the length of the slidign window with zeros at both ends
  smooth2=filter(temp2,hw2) # smooth the temp vector using the hann window
  smooth2=smooth2[(win2+1):(length(smooth2)-win2)] # trim smooth to correct length
  sae2.1[j+6]=smooth2
}

colnames(sae2.1)=c("Chr", "Pos", "Watson", "Crick", "Total", "watson.s", "crick.s", "total.s","watson.s2", "crick.s2", "total.s2")

##########################################################################################################################################################
# Plot boundaries:
#Adjust y lim here
plot(sae2.1$Pos,sae2.1$total.s/Mreads[k], type="n", xlim=c(xl1,xl2), ylim=c(0,1500), 
     ylab=paste(c(DSBListNames[k]))) #plot the start histogram


# Broad Overlays: Uncoment other lines if aiming to plot W and C separately
#lines(sae2.1$Pos,sae2.1$watson.s2*sae2scalar2/Mreads[k], type="h", xlim=c(xl1,xl2), col="lightcoral") #plot the start histogram
#lines(sae2.1$Pos,-sae2.1$crick.s2*sae2scalar2/Mreads[k], type="h", xlim=c(xl1,xl2), col="lightblue") #plot the start histogram
#lines(sae2.1$Pos,0.5*sae2.1$total.s2*sae2scalar2/Mreads[k]-1500, type="l", xlim=c(xl1,xl2), col="black") #plot the start histogram
lines(sae2.1$Pos,sae2.1$total.s2*sae2scalar2/Mreads[k], type="h", xlim=c(xl1,xl2), col=col.vec[k]) #plot the start histogram
}
##########################################################################################################################################################
#Now plot the gene datatrack
#First subset the relevant data
genes=AllElementsDUB #First make a copy of the ALLElements table
genes=subset(genes,chr==chrom & start>(xl1-10000) & stop<(xl2+10000)) #Make a sub-table of ALLElements where chr = 1 and has limits just beyond plot range
genes=subset(genes,type=="gene") #Make a sub-table of ALLElements

#Now perform the plot
plot(sae2.1$range,sae2.1$filtered, xaxt="n",yaxt="n",type="n", ylab=paste("Genes"),cex.lab=1.5,font=2, xlim=c(xl1,xl2), ylim=c(-100,120),axes=F) #set up empty plot
text((xl1+xl2)/2,-80, labels=paste("Chromosome",chrom, "/",orf,"/ Range",xl1,"to",xl2,"bp / Hann", win, "/ Y-Scalar",scalar, "/ backgroundsubtract", bgsubtract), cex.lab=1.4)
#abline(v=CEN[chr], col="dimgrey") 


##########################################################################################################################################################
########### TO PLOT GENE NAMES #############
##########################################################################################################################################################
# Following module draws arrows for each element
xrange=xl2-xl1
ahead=xrange/25 #make arrowhead length proportional to plot range
ahead[(ahead>500)]=500 #limit max length to 500
av=75 #arrow vertical location relative to plot dimensions
ahw=15 #arrow/head width
genesW=subset(genes,genename=="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
for (i in 1:nrow(genesW)){
  polygon(c(genesW[i,"start"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead,genesW[i,"stop"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead, genesW[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="wheat", border="wheat4", lty=2)
  text((genesW[i,"start"]+genesW[i,"stop"])/2,av, font=3, genesW[i,"sysname"], cex=0.9) }
genesW=subset(genes,genename !="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
for (i in 1:nrow(genesW)){
  polygon(c(genesW[i,"start"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead,genesW[i,"stop"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead, genesW[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="wheat", border="wheat4")
  text((genesW[i,"start"]+genesW[i,"stop"])/2,av, font=3, genesW[i,"genename"], cex=0.9) }

av=25 #arrow vertical location for Crick genes relative to plot dimensions
genesC=subset(genes,genename=="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
for (i in 1:nrow(genesC)){
  polygon(c(genesC[i,"stop"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead,genesC[i,"start"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead, genesC[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="thistle", border ="thistle4", lty=2)
  text((genesC[i,"start"]+genesC[i,"stop"])/2,av, font=3, genesC[i,"sysname"], cex=0.9) }
genesC=subset(genes,genename !="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
for (i in 1:nrow(genesC)){
  polygon(c(genesC[i,"stop"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead,genesC[i,"start"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead, genesC[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="thistle", border ="thistle4")
  text((genesC[i,"start"]+genesC[i,"stop"])/2,av, font=3, genesC[i,"genename"], cex=0.9) }
##########################################################################################################################################################

