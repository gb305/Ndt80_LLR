##############################################################################
########################## HOTSPOT IDENTIFICATION ############################
##############################################################################
#Packages to load
.libPaths("/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(data.table) ; library(e1071) ; library(GenomicFeatures) ; library(VennDiagram)

#Load chromosome number and size
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES")
chromlengths <- read.table("yeast_chrom_sizes.txt")

#Set working directory to where the .rds file is (.rds containing fullmap info of all the strains used to identify hotspota)
setwd("/Users/ll381/Dropbox/TESIS_DOCTORAL/R/RFILES/TEMPLATE_STUDY/CUTOOF_0.193/tel1D_control_rec8D_rec8Dtel1D") # Set wd
parent.directory <- getwd()
files <- list.files(pattern = ".rds");files
file.to.load <- 2L
load(files[file.to.load])
exp.name <- substr(files[file.to.load], 0, nchar(files[file.to.load])-4) # name of rds file
samples <- length(DSBList) #Number of libraries where hotspots will be identified
Mreads <- sapply(DSBList, function(x) { sum(x$Watson + x$Crick)/ 1000000}) #Million reads of each library

##############################################################################
#Hotspots characteristic
  win = 201  ##hann window 
  hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
  scalar=1/sum(hw) #scalar [1]
  cutoff = 0.193  #Cutoff to threshold libraries ; 0.193 HpM used in Pan and Neeman ; 0.125 used for spo11-yf 
  
  #############
  output.list2 <- rep(list(NULL), samples)
  for (k in 1:samples){
    output.list <- rep(list(NULL),16)
  for (chrom in 1:16) {
    sae2.0=subset(DSBList[[k]], Chr==chrom)
    sae2.1 <- data.frame(Chr=chrom, Pos=(1:chromlengths[chrom,2])) # Creates expanded empty dataframe with Chr and Pos locations
    sae2.1 <- merge(sae2.1,sae2.0, all=TRUE) # Merge expanded empty dataframe with compressed sae2.1 dataframe
    sae2.1[is.na(sae2.1)] <- 0 # Convert all NA values to zero
    sae2.1$Total=sae2.1$Watson+sae2.1$Crick
    sae2.1$HpM <- sae2.1$Total/ Mreads[k] # convert Total to HpM
    sae2.1$Smooth = scalar*as.integer(filter(sae2.1$HpM, hw)) ##Smooth data
    sae2.1[is.na(sae2.1)] = 0L
    sae2.1 <- sae2.1[sae2.1$Smooth >= cutoff,]  ### HERE IS THE threshold filter
    gr <- GRanges(seqnames = sae2.1$Chr, ranges = IRanges(start = sae2.1$Pos, end = sae2.1$Pos))
    gr2 <- reduce(gr, ignore.strand = T)  ### REDUCE FUCTION CREATES ONE SEGMENT IN OVERLAPED REGIONS
    ##This merge + 100bp adjacent to the hotspot segment
    start(gr2) <- start(gr2) - 100
    end(gr2) <- end(gr2) + 100
    ##Reduce again overlapped segments in one
    gr2 <- reduce(gr2, ignore.strand = T)
    start(gr2) <- start(gr2) + 100
    end(gr2) <- end(gr2) - 100
    
    Chr <- as.numeric(as.character(seqnames(gr2)))
    Start <- as.numeric(start(gr2))
    End <- as.numeric(end(gr2))
    mask <- cbind.data.frame(Chr,Start,End)
    mask <- mask[(mask$End - mask$Start) >= 25,]  ## APPLY A LENGTH FILTER, NO HT IF IT IS < 25 bp
    gr2 <- GRanges(seqnames = mask$Chr, ranges = IRanges(start = mask$Start, end = mask$End)) 
   
    sae2.2 <- sae2.1[sae2.1$Total > 0,]
    gr3 <- GRanges(seqnames = sae2.2$Chr, ranges = IRanges(start = sae2.2$Pos, end = sae2.2$Pos))
    overlaps <- findOverlaps(gr2,gr3)
    df <- cbind.data.frame(Smooth = queryHits(overlaps), Nuc = subjectHits(overlaps))
  
   output <- data.frame(matrix(NA, ncol = 3, nrow = 10000)); 
   colnames(output) <- c("Chr", "Start", "End") #Creates a regular format Chr, Start and End column
   output$Chr <- chrom
  
  ###Filter HPs that dont have more than 25 reads
  for(i in 1:max(df$Smooth)) {
    temp <- df[df$Smooth == i,]
    if (nrow(temp) < 25) next
    output[i,"Start"] = start(gr3)[temp[1,"Nuc"]]
    output[i,"End"] = start(gr3)[temp[nrow(temp), "Nuc"]]
  }
  output <- output[!is.na(output$Start),]
    output.list[[chrom]] <- output
  }
    output.list <- do.call("rbind", output.list)
    output.list2[[k]] <- output.list
  }

  # This part combines the hotspots from all strains to give a single "mask" of all hotspots found in every strain
  output.bind <- do.call("rbind", output.list2)
  #Merge the hotspot that overlap X (We chose 12 since the minimum hotspot lenght defined was 25bp)
  gr4 <- GRanges(seqnames = output.bind$Chr, ranges = IRanges(start = output.bind$Start, end = output.bind$End))
  start(gr4) <- start(gr4) + 12
  end(gr4) <- end(gr4) - 12
  gr4 <- reduce(gr4, ignore.strand = T)  ### REDUCE FUCTION CREATES ONE SEGMENT IN OVERLAPED REGIONS
  #Correct again
  start(gr4) <- start(gr4) - 12
  end(gr4) <- end(gr4) + 12

  Chr <- as.numeric(as.character(seqnames(gr4)))
  Start <- as.numeric(start(gr4))
  End <- as.numeric(end(gr4))
  HS.store <- cbind.data.frame(Chr,Start,End)

  ############################################################   
  ## Save new template
  write.table(HS.store, file = paste0(cutoff, "all_template.txt"), row.names = F, col.names = T, quote = F, sep = "\t")  
  save(HS.store, file = "all_template.rds")  
  ##Save the data in a .rds file
  save(output.list2, file = paste0(cutoff,"HpM_all_template.rds"))
  