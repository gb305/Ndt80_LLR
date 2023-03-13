#FUNCTIONS

CCfoldHS.calc<-function(p.dir=getwd(),s.dir,c.dir,coln="NormHpMChr",readfilename="Hotspot.Table",constant=0){


  setwd(c.dir);controls <- list.files(pattern = readfilename);print(controls)

  setwd(s.dir);files <- list.files(pattern = readfilename);print(files)
  setwd(p.dir)



  hs.renamer<-function(input){
    output=strsplit(input, "Hotspot.Table.")[[1]][2]
    output = strsplit(output, "[, _ -]+")[[1]][1]
    output = strsplit(output, ".txt")[[1]][1]
    print(output)
    return(output)
  }
  packages.check("data.table")
  for (c in 1:length(controls)){
    for (i in 1:length(files)){
      setwd(s.dir);sv=fread(files[i])
      setwd(c.dir);cv=fread(controls[c])
      if (files[i]!=controls[c]){
        sv$log2fold=log2((sv[[coln]]+constant)/(cv[[coln]]+constant))
        fg=data.frame(Chr=sv$Chr,Midpoint=sv$Midpoint,log2fold=sv$log2fold,minname=hs.renamer(controls[c]),plusname=hs.renamer(files[i]))
        print(paste0(hs.renamer(files[i]),"/",hs.renamer(controls[c])))
        fg[sapply(fg, is.infinite)] <- NA
        setwd(p.dir);write.table(fg,file=paste0(hs.renamer(files[i]),"-",hs.renamer(controls[c]),"_",coln,".txt"),quote = FALSE, sep="\t",row.names = FALSE)
      }
    }
  }
}
CCfoldHS.smoother<-function(input,sc=80,si=100,txtout=F){
  chrSize=c(230218,813184,320870,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
  outputdf=fread(input)
 outputdf<-na.omit(outputdf)
  chrpl=list()
  for (c in 1:16){
    print(c)
    chrp=subset(outputdf, Chr ==c)

    chrp$log2fold[chrp$log2fold == -Inf] <- 0
    chrp$log2fold[chrp$log2fold == Inf] <- 0
    loessMod10 <- loess(control=loess.control(surface="direct"),formula= as.numeric(log2fold) ~ as.numeric(Midpoint), data=chrp, span=sc/nrow(chrp))
    chrpo <- data.frame(Midpoint = seq(1,chrSize[c],si)) #change to 100bp
    chrpl[[c]]=transform(chrpo, log2fold = predict(loessMod10, chrpo))# 10% smoothing span
    chrpl[[c]]$Chr=c

  }
  out.smooth<- do.call("rbind",  chrpl)
  out.smooth$log2fold=round(out.smooth$log2fold,3)
  dir.create("smoothed")

  # write.table(jk,file=paste0("smoothed/smooth_",input), sep = "\t", row.names = F, quote = F)
  save(out.smooth,file=paste0("smoothed/smooth_",substr(input, 1, nchar(input)-4) ,"_sc",sc,".Rbin"), compress=T)
  if(txtout==T){write.table(out.smooth, file = paste0("smoothed/smooth_",substr(input, 1, nchar(input)-4) ,"_sc",sc,".txt"),  sep = "\t", row.names = F, quote = F)}
}

CCfoldHS.binning <-function(files,i=1,binwidth.v=c(10,25,50,100),work.d=getwd(),si=100,plotmode=T,ARSline=F,ylims=2,control.col="lightblue",mutant.col="red",usesmooth=T) {
  if (usesmooth==T){HS <- loadRData(files[i])}else{HS=fread(files[i]);plotmode=F}


      sfiles = strsplit(files[i], ".txt")[[1]][1]
      if (usesmooth==T){sfiles = strsplit(sfiles, "[, _ ]+")[[1]][2]}
  s1=strsplit(sfiles, "[, _ -]+")[[1]][1]
  s2=strsplit(sfiles, "[, _ -]+")[[1]][2]
  for (v in 1:length(binwidth.v)){
    binwidth=binwidth.v[v]
    # Mreads[k]<-sum(temp$Watson + temp$Crick)/ 100000

    convert.dsblist <- function(input){
      output <- data.frame("Chr" = seqnames(input), "start" = start(input), "end" = end(input), "log2fold" = input$log2fold)    #convert back to DSBList format do num of reads in each stacked bar chart.
      return (output)
    }

    # genome <- BSgenome.Scerevisiae.UCSC.sacCer3
    # Chr=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
    # Chr=as.character(c(1:16))
    chrSize=c("chrI"=230138,"chrII"=813197,"chrIII"=316643,"chrIV"=1531994,"chrV"=439575,"chrVI"=270161,"chrVII"=1090940,"chrVIII"=562643,"chrIX"=439888,"chrX"=745751,"chrXI"=666816,"chrXII"=1078177,"chrXIII"=924431,"chrXIV"=784333,"chrXV"=1091291,"chrXVI"=948066)
    chrSizeChr=c("1"=230218,"2"=813184,"3"=320870,"4"=1531933,"5"=576874,"6"=270161,"7"=1090940,"8"=562643,"9"=439888,"10"=745751,"11"=666816,"12"=1078177,"13"=924431,"14"=784333,"15"=1091291,"16"=948066)

    chrSize=c(230218,813184,320870,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
    hj=list()
    binwidthbp=binwidth*1000
    for (c in 1:16){
      hj[[c]]=subset(HS,Chr==c)

      # hj[[c]]=tail(hj[[c]],-round(chrSize[c]%%binwidthbp/(si*2)))

      # hj[[c]]=head(hj[[c]],-round(chrSize[c]%%binwidthbp/(si*2)))
       chrSizeChr[c]=chrSize[c] %/% binwidthbp*binwidthbp
    }
    print(chrSizeChr[c])
    cd=NULL
    cd=do.call(rbind, hj)
    cd$Chr=as.character(cd$Chr)
    if (usesmooth==T){cd$log2fold=cd$log2fold*si}else{cd$log2fold=cd$log2fold*(sum(chrSize)/nrow(cd))}
    cd$log2fold[cd$log2fold == -Inf] <- 0
    cd$log2fold[cd$log2fold == Inf] <- 0
    cd[is.na(cd)] <- 0
    CCseq_gr <- makeGRangesFromDataFrame(cd,                           # generate a GRange for your CCseq data
                                         keep.extra.columns=T,
                                         ignore.strand=T,
                                         seqinfo=NULL,
                                         seqnames.field=c("Chr"),
                                         start.field="Midpoint",
                                         end.field=c("Midpoint"),
                                         starts.in.df.are.0based=FALSE)

    cd=NULL
    gr.data.cov <- GenomicRanges::coverage(CCseq_gr,weight = "log2fold")
    gr.windows=tileGenome(chrSizeChr,tilewidth=binwidth*1000,cut.last.tile.in.chrom=T)
    # ranges(gr.windows)=IRanges(start=230218%%100000/2:230218,width=99999)
    # gr <- GRanges(seqnames = c(1:2), ranges = IRanges(start = 1:200000, width = 1000))
    # startx=seq(1,100001,100000)
    # x <- Seqinfo(seqnames=as.character(c(1:16)),
    # seqlengths=chrSize,
    # genome="toy")
    # gr1 <- GRanges("3:15-25", seqinfo=x)
    # gr1 <- GRanges("3:15-25", seqinfo=x)
    # gr <- GRanges(seqnames = 1, ranges = IRanges(start = 1:20000, Endwidth = 2000))
    # gr <- GRanges(seqnames = c("chr1", "chrU345"),
    # ranges = IRanges(start = startx, end = startx+100000))
    # gr.windows=tileGenome(seqinfo(genome),tilewidth=binwidth*1000,cut.last.tile.in.chrom=T)
    # gr.windows=unlist(gr.windows)
    seqlevels(gr.windows, pruning.mode="coarse") <- names(gr.data.cov)

    gr.data.binnedAvg <- binnedAverage(gr.windows, gr.data.cov, "log2fold")
    df=convert.dsblist(gr.data.binnedAvg)
    gr.data.cov=NULL
    CCseq_gr=NULL
    gc()
    print(df[1,4])

    df$mid=(df$end-df$start)/2+df$start
    temp=c()
    temp= ifelse(df$log2fold >0,paste0(s1),paste0(s2))
    temp2= ifelse(df$log2fold ==0,NA,df$log2fold)
    temp3= ifelse(is.na(temp2),temp2,temp)
    df$Library=temp3
    # df$Chr=as.numeric(df$Chr)
    dir.create("binData");dir.create(paste0("binData/",binwidth,"Kb"))
    dir.create("binPlots");dir.create(paste0("binPlots/",binwidth,"Kb"))
    write.table(df,file=paste0("binData/",binwidth,"Kb/",binwidth,"Kb_",substr(files[i], 1, nchar(files[i])-5),".txt"), sep = "\t", row.names = F, quote = F)
    if(plotmode==T){
    centro=CCAnnotate(Cer3H4L2_centro)
    CEN=round((centro$chromEnd+centro$chromStart)/2)
    # CEN=c(151523, 238265, 114443, 449766,152045, 148568, 496979,105644, 355687, 436366,440187, 150887, 268090, 628816, 326643, 556015) # CENtromere positions
    c=10
    options(scipen=10000)
    binplot=df
    print(colnames(binplot))
    # dir.create(paste0("binPlots/",binwidth))
    multi.page=list()
    ARS=CCAnnotate(ARS_consensus)

    test <- CCAnnotate(Cer3H4L2_AllElementsDUB)
    levels(test$type)
    ARS <- test[test$type == "ARS_consensus_sequence",]
    parent.directory <- dirname(work.d)
    if (usesmooth==T){setwd(parent.directory)}
    unsmoothed=list.files(pattern=sfiles)
    library(data.table)
    df=fread(unsmoothed[1])
    setwd(work.d)
    # test <- test[test$chr == 11,]
    for (c in 1:16){
      binchr=subset(binplot,Chr==c)
      print(colnames(binchr))
      toplot=subset(df,Chr==c)
      smooth=subset(HS,Chr==c)
      ARSchr=subset(ARS,chr==c)
      arsmidpoint=as.vector((ARSchr$stop-ARSchr$start)+ARSchr$start)
      ARSv=as.vector(ARSchr$Pos)
      ARSv=arsmidpoint

      # toplot=head(toplot,-1)
      # toplot$mid=toplot$mid/1000
      # toplot$start=toplot$start/1000
      # toplot$end=toplot$end/1000
      if (ARSline==F){ARSv=-1000}
      toplot$Midpoint=toplot$Midpoint/1000
      smooth$Midpoint=smooth$Midpoint/1000
      binchr$mid=binchr$mid/1000
      multi.page[[c]]=ggplot(smooth,aes(Midpoint,log2fold))+
        geom_line()+
        geom_bar(data=binchr,aes(mid,log2fold,fill=Library),stat="identity",width=binwidth,alpha=0.3)+
         geom_vline(xintercept = ARSv/1000)+

        annotate(geom = "label",x=CEN[c]/1000,y=-3, label = "Centromere", show.legend = FALSE,alpha=0.2)+
        # annotate(geom = "label",x=ARSv/1000,y=2, label = "ARS", show.legend = FALSE,alpha=0.2)+

         geom_histogram(stat="identity",data=toplot, aes(Midpoint,log2fold),width=0.5)+

        # geom_line(data=smooth, aes(Midpoint,log2fold))+
        # geom_vline(xintercept = ARSv)+ geom_line(HS[[v]])
        # geom_vline(xintercept = CEN[c])+
        theme_bw(base_size = 18) +
        # +theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        ylim(-ylims,ylims)+xlim(0,max(smooth$Midpoint))+

        labs(x = paste0("Chromosome ", c, " (Kb)"),y=paste0("Fold log2(",strainTimeSplit(s1),"/",strainTimeSplit(s2),")")) + scale_fill_manual(values = c(control.col, mutant.col))

      # +scale_x_continuous(breaks=seq(0,max(toplot$end),100))

    }
    ggexport(multi.page, filename = (paste0("binPlots/",binwidth,"Kb/",binwidth,"Kb_", substr(files[i], 1, nchar(files[i])-5) ,".pdf")),width = 14,height = 7)
    # chrSize=c("1"=230218,"2"=813184,"3"=320870,"4"=1531933,"5"=576874,"6"=270161,"7"=1090940,"8"=562643,"9"=439888,"10"=745751,"11"=666816,"12"=1078177,"13"=924431,"14"=784333,"15"=1091291,"16"=948066)
    }
  }
}
CCfoldHS.binpileup <-function(p.dir=getwd(),binwidth.v=c(100,50,25,10),log2FoldMax=2,centromere=T,EXT=".png",invert.col=F,fixmax=F,log2foldoffset=0,colhigh="red",collow="navy",usesmooth=T){
  packages.check("ggplot2");packages.check("ggpubr");packages.check("lemon")
  for (i in 1:length(binwidth.v)){
    setwd(paste0(p.dir,"/binData/",binwidth.v[i],"Kb"));d.dir=getwd()
    readfilename=paste0(binwidth.v[i],"Kb")
    files = list.files(pattern=readfilename)

    for (f in 1:length(files)){
      sfiles = strsplit(files[f], ".txt")[[1]][1]
      if (usesmooth==F){ti=1}else{ti=0}
      s1=strsplit(sfiles, "[, _ -]+")[[1]][3-ti]
      s2=strsplit(sfiles, "[, _ -]+")[[1]][4-ti]
      coln=strsplit(sfiles, "[, _ -]+")[[1]][5-ti]
      HS <- fread(paste0(files[f]))
      HS$log2fold <- HS$log2fold - log2foldoffset
      if (fixmax==T){
        HS$log2fold=ifelse(HS$log2fold>=log2FoldMax,log2FoldMax,HS$log2fold)
        HS$log2fold=ifelse(HS$log2fold<=-log2FoldMax,-log2FoldMax,HS$log2fold)

      }

      if (invert.col==T){
        temptg=s1
        s1=s2
        s2=temptg
        HS$log2fold=-HS$log2fold
      }
      binwidth=binwidth.v[i]

      HS$length=HS$end-HS$start
      hsl=list()

      CEN=c(151523, 238265, 114443, 449766,152045, 148568, 496979,105644, 355687, 436366,440187, 150887, 268090, 628816, 326643, 556015) # CENtromere positions
      if (centromere==F){CEN=rep(1,16)}
      for (c in 1:16){
        HSt=subset(HS, Chr ==c)
        e=subset(HSt,end<CEN[c]+binwidth*1000 &start>CEN[c]-binwidth*1000)
        cenr=round(e$mid)
        HSt$cen=round((HSt$mid-cenr)/1000)
        hsl[[c]]=HSt
      }
      HSc=do.call("rbind", hsl)

      max=max(HSc$end)
      HSc$cons=1


      chrSize=c("1"=230218,"2"=813184,"3"=320870,"4"=1531933,"5"=576874,"6"=270161,"7"=1090940,"8"=562643,"9"=439888,"10"=745751,"11"=666816,"12"=1078177,"13"=924431,"14"=784333,"15"=1091291,"16"=948066)
      chrsorder =c(4, 15,7,12,16,13,2,14,10,11,5,8,9,3,6,1)
      options(scipen=10000)

      xlims.max=1200
      xlims.min=-xlims.max
      axis.labal.max=1100
      axis.labal.min=-600
      x.label="Position relative to centromere (kb)"
      titlefront="CEN"
      if (centromere==F){
        xlims.max=1800
        xlims.min=0
        axis.labal.max=1500
        axis.labal.min=0
        x.label="Position relative to the left telomere (kb)"
        titlefront="TEL"

      }

      HSc=na.omit(HSc)
      HSc$Chr <- factor(HSc$Chr,      # Reordering group factor levels
                        levels = rev(chrsorder))
      g=ggplot(HSc,aes(x=cen,y=cons,fill=log2fold))+geom_col(width=binwidth)+facet_grid(margins=F,rows=vars(Chr),switch="y")+xlim(xlims.min,xlims.max)+scale_fill_gradient2(
        low = collow, high = colhigh,limits=c(-log2FoldMax,log2FoldMax),breaks=c(-log2FoldMax,0, log2FoldMax),labels=c(paste0(-log2FoldMax,"  ",strainTimeSplit(s2)),0,paste0(log2FoldMax,"  ",strainTimeSplit(s1))))+
        theme_classic2(base_size = 15)+theme(panel.background=element_rect("lightgrey"),strip.text.y.left = element_text(angle = 0),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.line.y=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(),plot.title = element_text(size=12))+
        xlab(x.label)+labs(title=paste0(coln," Log2 fold change (",strainTimeSplit(s1),"/",strainTimeSplit(s2),") binned at ",binwidth,"kb"))+
        coord_capped_cart(xlim=c(axis.labal.min,axis.labal.max),bottom='both')



      setwd(p.dir)
      dir.create(paste0("binPlots/",binwidth.v[i],"Kb/",titlefront,"pileup"))
      ggsave(g, filename = (paste0("binPlots/",binwidth.v[i],"Kb/",titlefront,"pileup/",titlefront,"pileup_",binwidth,"Kb_",s1,"-",s2,"_",coln,EXT)),width = 8,height = 6)
      # ggexport(g, filename = (paste0("binPlots/",binwidth.v[i],"Kb/",titlefront,"pileup_",binwidth,"Kb_",s1,"-",s2,"_",coln,".pdf")),width = 8,height = 6)

      setwd(d.dir)
    }
  }
  setwd(p.dir)
}

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
 
