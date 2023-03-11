##############################HOTSPOT_ANALYSIS###############################
library(data.table) ; library(e1071) ; library(GenomicFeatures) ; library(VennDiagram)
##############################################################################
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES")
chromlengths <- read.table("yeast_chrom_sizes.txt")
##############################################################################

#Hotspot template comparison

setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES")
Pan <- read.delim("Pan.Hotspots.IGR.SacCer3_H4L2_2016.08.10a.txt")
Neeman <- read.delim("Neeman_template.txt")
Neale <- read.delim("0.193V2_template_UPDATE_3473HS.txt")
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/CUTOOF_0.193/tel1D_control/spo11-yf")
Spo11_yf <- read.delim("0.125_Hotspot.Table.spo11_yf_NorDNA.txt")

#template1  = first template (Normally Neale template); template2 = the other template to compare 
template1 =GRanges(seqnames = Neale$Chr, ranges = IRanges(start = Neale$Start, end = Neale$End))
template2 <- GRanges(seqnames = Neeman$Chr, ranges = IRanges(start = Neeman$Start, end = Neeman$End)) # Comment if not used
template2 <- GRanges(seqnames = Pan$CHROM, ranges = IRanges(start = Pan$HS_START, end = Pan$HS_END)) # Comment if not used
template2 <- GRanges(seqnames = Spo11_yf$Chr, ranges = IRanges(start = Spo11_yf$Start, end = Spo11_yf$End)) # Comment if not used

# Find number of hotspot that overlap in position 
overlap.pan <- findOverlaps(template1, template2)
Pan2 <- length(template2)
Control <- length(template1)
Shared2 <- length(overlap.pan)

#Plot Venn diagram
plot.new()
draw.pairwise.venn(Control, Pan2, Shared2, category = c(" "," "),
                   cat.pos = c(10,100), cex=1.5,cat.cex = 1.5,
                   fontface = "bold",
                   col = c("lightblue4","ivory4"), 
                   fill = c("lightblue2", "ivory3"),
                   fontfamily = c("Helvetica") )

Pan2 ;Control ;Shared2

#COLOURS USED IN THIS STUDY:
#FOR NEALE AND NEEMAN line "lightblue4","peachpuff3" ; lightblue2", "peachpuff" filled
#FOR NEALE AND PAN  col = c("lightblue4","lightpink4"), fill = c("lightblue2", "lightpink4"),
#FOR spo11-yf col = c("lightblue4","ivory4"),  fill = c("lightblue2", "ivory3")

##OVERLAPPING PERCENTAGE
NotOVerlaped <- (Pan2-Shared2) + (Control-Shared2)
total <- NotOVerlaped+Shared2
Percentaje <- Shared2/total*100
NotOVerlaped ; total ; Percentaje #Percentage of hotspot overlap

###############################################################
####### Frequency and strength of exclusive hotspots 
#thius section finds out the frequency and strength of the exclusive hotspots present between any two templates
#For neale template:
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/CUTOOF_0.193/tel1D_control")
control = read.delim("V2T_Hotspot.Table.sae2Dndt80D.txt")
# #For neeman template
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/CUTOOF_0.193/tel1D_control/NEEMAN")
samp <- read.delim("Hotspot.Table.Cer3H4L2_NeemanTel1D_45566hAverage_-DC>32bp.txt") #Hotspot table
# ##samp <- read.delim("Neeman_template.txt")  #Hotspot template

#For Pan
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/TEMPLATE_ANALYISIS/PAN")
samp <- read.delim("TEMPLATEPAN_Hotspot.Table.pan.txt") #Hotspot table
#samp <- read.delim("Pan.Hotspots.IGR.SacCer3_H4L2_2016.08.10a.txt") #Hotspot template

#SPO11-yf
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/CUTOOF_0.193/tel1D_control/spo11-yf")
samp <- read.delim("0.125_Hotspot.Table.spo11_yf_NorDNA.txt")

template1.gr <- GRanges(seqnames = control$Chr, ranges = IRanges(start = control$Start, end = control$End), mcols = control[,4:18])
samp.gr <- GRanges(seqnames = samp$Chr, ranges = IRanges(start = samp$Start, end = samp$End), mcols = samp[,4:ncol(samp)]) #FOR THE OTHER DATA

#query #subject
overlap_n2 <- findOverlaps(template1.gr, samp.gr)

# BEWARE! If you use the numbers below tomake a venn diagram, it will not be accurate, because the two sets do not BIJECT
Shared2 <- length(overlap_n2)
samp2 <- length(samp.gr)
Control <- length(template1.gr)
samp2 ;Control ;Shared2

#these are the the real hotspotds that are exclusive tothe control. 
# Note the number is different to the number in the Venn diagram,because the venn diagramis inaccurate.
# This is because there are occasions where region defined as ONE hotspot in one template, is defined as TWO (or more) in another. 
# This means that there is not one number of shared hotspots. 
# Put another way, the number of control hotspots that is shared with samp is not the same as the number of samp hotspots that are shared with control.
# I.E there is not a ONE-TO-ONE correspondence between the sets.

onlyincontrol <- control[-queryHits(overlap_n2),] #FIND HOTSPOTS ONLY PRESENT IN CONTROL
onlyinsample <- samp[-subjectHits(overlap_n2),] #FIND HOTSPOTS ONLY PRESENT IN SAMPLE
controlshared <- control[queryHits(overlap_n2),] #FIND HOTSPOTS SHARED IN BOTH CONTROL AND SAMPLE

# Plot the subset of exclusive hotspots in Template 1
plot1 <- density(log2(control$NormHpM))
plot2 <- density(log2(onlyincontrol$NormHpM))
plot(plot1$x, plot1$y*plot1$n, type="l", xlab = "Log2(NormHpM)", ylab = "Frequency", col = "lightblue4",lwd=2,
     ylim= c(0,1000),  xlim= c(0,15)) # the dsitribution of hotspot strength in the whole control set
lines(plot2$x, plot2$y*plot2$n, type="l", col = "lightblue3",lwd=2) # the distribution of hotspot strength in the subset that is EXCLUSIVE to control (vs sample)

# Plot the subset of exclusive hotspots in Template 2
plot3 <- density(log2(samp$TotalHpM))
plot4 <- density(log2(onlyinsample$TotalHpM))
plot(plot3$x, plot3$y*plot3$n, type="l", xlab = "Log2(TotalHpM)", col = "lightpink4",
     ylab = "Frequency", ylim= c(0,1000),xlim= c(0,15),
     lwd=2) # the dsitribution of hotspot strength in the whole sample set
lines(plot4$x, plot4$y*plot4$n, type="l", col = "lightpink3",lwd=2) # the distribution of hotspot strength in the subset that is EXCLUSIVE to sample (vs control)

