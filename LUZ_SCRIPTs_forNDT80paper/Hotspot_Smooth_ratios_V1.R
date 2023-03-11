###########################################################################
##########################  HOTSPOT_ RATIOS  ########################
###########################################################################
#Load chromosome number and length
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES")
chromlengths <- read.table("yeast_chrom_sizes.txt")
chromlengths$kb <- round(chromlengths[,2]/1000, digits=1)

#Setwt to Hotspot.table files
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/CUTOOF_0.193/tel1D_control")
parent.directory <- getwd()
setwd(parent.directory); 

# Centromere positions
CEN=c(151523.5, 238265, 114443, 449766,152045.5, 148568.5, 496979,105644.5,
      355687, 436366,440187.5, 150887.5, 268090, 628816.5, 326643, 556015) 

##Load data - set each one in a separate variable
files = list.files(pattern="Hotspot.Table");files
DSBListNames = substr(files,21, nchar(files)-4) ; DSBListNames #Assign names
#Load files in the desired order
c1 <- read.table(files[1], sep = "\t", header=TRUE)
c2 <- read.table(files[2], sep = "\t", header=TRUE)
c3 <- read.table(files[3], sep = "\t", header=TRUE)
c4 <- read.table(files[4], sep = "\t", header=TRUE)

contitions <- c("c1","c2","c3","c4")
names = DSBListNames
names = c(DSBListNames[1],DSBListNames[2],DSBListNames[3],DSBListNames[4])
names(contitions) <- names ; names2 <- names ; contitions 

##################################################################################################################################
##Add a constant value  to NormHpM in order to avoid big differences from weak hotspots
x <- 5 #Adjust constant value
c1$NormHpX_PLUS<-c1$NormHpM + x ; c2$NormHpX_PLUS<-c2$NormHpM + x ; 
c3$NormHpX_PLUS<-c3$NormHpM + x ; c4$NormHpX_PLUS<-c4$NormHpM + x

##Generate a new NormHpMChr column in each strain
rowA=0 ; rowB=0; rowC=0 ; rowD=0
#For control
for (i in 1:16){
  aChr=subset(c1,Chr==i)
  aChrDSBs=nrow(aChr) # Number of hotspots per chromosome (rows)
  aChrTotal=sum(aChr[1:aChrDSBs,"NormHpX_PLUS"])
  c1[(rowA+1):(aChrDSBs+rowA),"NormHpChr_PLUS"]=aChr[1:aChrDSBs,"NormHpX_PLUS"]/aChrTotal*1000000
  rowA=rowA+aChrDSBs
}
#For rec8

for (z in 1:16){
  aChr=subset(c2,Chr==z)
  aChrDSBs=nrow(aChr) # Number of hotspots per chromosome (rows)
  aChrTotal=sum(aChr[1:aChrDSBs,"NormHpX_PLUS"])
  c2[(rowB+1):(aChrDSBs+rowB),"NormHpChr_PLUS"]=aChr[1:aChrDSBs,"NormHpX_PLUS"]/aChrTotal*1000000
  rowB=rowB+aChrDSBs
}
#For tel1

for (c in 1:16){
  aChr=subset(c3,Chr==c)
  aChrDSBs=nrow(aChr) # Number of hotspots per chromosome (rows)
  aChrTotal=sum(aChr[1:aChrDSBs,"NormHpX_PLUS"])
  c3[(rowC+1):(aChrDSBs+rowC),"NormHpChr_PLUS"]=aChr[1:aChrDSBs,"NormHpX_PLUS"]/aChrTotal*1000000
  rowC=rowC+aChrDSBs
}
#For tel1rec8

for (v in 1:16){
  aChr=subset(c4,Chr==v)
  aChrDSBs=nrow(aChr) # Number of hotspots per chromosome (rows)
  aChrTotal=sum(aChr[1:aChrDSBs,"NormHpX_PLUS"])
  c4[(rowD+1):(aChrDSBs+rowD),"NormHpChr_PLUS"]=aChr[1:aChrDSBs,"NormHpX_PLUS"]/aChrTotal*1000000
  rowD=rowD+aChrDSBs
}
 

##Generate a large list that contains all the samples together
data <- list(c1,c2,c3,c4); names(data) <- names  

##Samples positions in the data list
sample1 <- 1 ; sample2 <- 2 ;
sample3 <- 3 ; sample4 <- 4
#####SPECIFY SAMPLE NUMBERS TO COMPARE
contitions 
condition1 <- sample1 #####DIVIDEND 
condition2 <- sample2 #####DIVISOR

#########################################################################################################
##########RUN FROM HERE

##Folder creation
##Set first working directory
setwd(parent.directory); 
ifelse(!dir.exists(file.path(parent.directory, "RATIOS_Chr")), dir.create(file.path(parent.directory, "RATIOS_Chr")), FALSE); 
setwd(file.path(parent.directory, "RATIOS_Chr")); 
parent.directory2 <- getwd()
setwd(parent.directory2); 
ifelse(!dir.exists(file.path(parent.directory2,paste0("Norm_HpM_",x))), dir.create(file.path(parent.directory2,paste0("Norm_HpM_",x))), FALSE); 
setwd(file.path(parent.directory2, paste0("Norm_HpM_",x))); 
directory3 <- getwd()
dir.create(file.path(directory3, paste0(names[condition1],"_vs_",names[condition2],"_PLUS_",x)), showWarnings = FALSE); 
setwd(file.path(directory3, paste0(names[condition1],"_vs_",names[condition2],"_PLUS_",x)))
directory4 <- getwd()

###RATIO - Specify if NormHpMChr or NormHpM 
data[[condition1]]$Log2Ratio_300 <- log2(data[[condition1]]$NormHpX_PLUS/data[[condition2]]$NormHpX_PLUS)
data[[condition1]]$Ratio_300 <- data[[condition1]]$NormHpX_PLUS/data[[condition2]]$NormHpX_PLUS

##Creation of a color vector, same samples lenght and position
col.vec <- c("dodgerblue4","dodgerblue1","goldenrod1","goldenrod4") #±Ndt80
#col.vec <- c("dodgerblue3","goldenrod3","dodgerblue3","goldenrod3") #NEEMAN strains
col.vec2 <- c(col.vec[condition2],col.vec[condition1]) #Assign colour to strain

pdf(file = paste0("Log2Ratio_NormHpM_",x,".pdf"), width = 6, height = 5)
for (b in 1:16) {
  chr <- b
  sample_chr <- data[[condition1]][data[[condition1]]$Chr == chr,]
  sample_chr <- na.omit(sample_chr)
  chr.length <- chromlengths[chr,2]
plot(1:chr.length, rep(0, chr.length),  type = "l", ylim = c(-2, 2), xlab = "Chromosome coordinates (bp)",
     ylab= "Log2 Ratio", main = paste0("Chromosome_", chr,"  ",chromlengths[b,3],"kb"))
lines(sample_chr$Midpoint, sample_chr$Log2Ratio_300, type = "h",lwd = 1,col = col.vec2[((sample_chr$Log2Ratio_300>0)+1)]) 
lines(smooth.spline(sample_chr$Midpoint, sample_chr$Log2Ratio_300,spar = 0.7), type = "l",col = "black", lwd = 4) 
points(CEN[chr], 0, col="black", lwd= 0.3, pch = 20)
}
dev.off()