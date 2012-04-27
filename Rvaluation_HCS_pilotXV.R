#requirements
require(plyr)
require(ggplot2)
require(gdata)

#trim transfection
kIn1$Transfection <- trim(kIn1$Transfection)

#exclude bad plates
kIn1.subset <- subset(kIn1, plate < 7 & plate > 0)

#calculate additional scores
kIn1.subset <- ddply(kIn1.subset, .(plate), transform, zscore_whole_plate=(CTblue - mean(CTblue)) / sd(CTblue),
                     rzscore_only_controls=zscore,
                     rzscore_whole_plate=(CTblue - median(CTblue)) / mad(CTblue), 
                     zscore_only_inner=my.zscore(plateColumn, CTblue), 
                     rzscore_only_inner=my.zscore(plateColumn, CTblue, robust=T))

#plot it all
my.createAllPlots <- function(kIn1.subset, margin=2, Bscore.margin=3, title="Pilot XV", method="mad", dir=NA, onlyLower=F, makeplot=T)
{
  #ask for target directory if not set
  if(is.na(dir) && makeplot)
  {
    dir <- choose.dir()
  }
  
  #iterate normalization types
  require(foreach)
  
  normalizationMethods <- c("CTblue", "Bscore", "poc", "npi",
                            "zscore_whole_plate", "rzscore_whole_plate",
                            "zscore_only_inner", "rzscore_only_inner")
  
  significantHitList <- data.frame()
  
  foreach(signalType = normalizationMethods) %do%
  {
    #set signal column
    kIn1.subset$signal <- kIn1.subset[[signalType]]
    
    #zscore-based?
    currentMethod <- method
    scoreBased <- F
    
    if(grepl("zscore", signalType))
    {
      scoreBased <- T
      currentMethod <- "none"
    }
    
    if(grepl("Bscore", signalType))
    {
      margin <- Bscore.margin
    }
    
    if(makeplot){
      #plot everything in one pdf file (per normalization)
      pdf(file=paste(dir, "/", signalType, ".pdf", sep=""), width=11.69, height=8.27, paper="a4r")
      
      my.heatmap(kIn1.subset, signalType, margin=margin, method=currentMethod)
      my.boxplot(kIn1.subset, signalType, plotCV=!scoreBased)
      my.plateComparison(kIn1.subset, signalType)
      my.well.to.well(kIn1.subset, signalType)
      
      #if zscore-based we can easily add a qqplot
      if(scoreBased)
      {
        my.qqplot(kIn1.subset, signalType, margin=margin)
      }
      dev.off()
    }
    
    #hitlists
    my.data <- kIn1.subset
    my.data.outliers <- cbind(subset(my.outliers(my.data, method=currentMethod, margin=margin, onlyLower=onlyLower), 
                               select=c(plate, plateRow, plateColumn, Transfection, signal)), normalization=signalType)
    
    if(dim(significantHitList)[1] == 0)
    {
      significantHitList <- my.data.outliers
    }
    else
    {
      significantHitList <- rbind(significantHitList, my.data.outliers)
    }
  }
  if(makeplot)
  {
    #compare normalizations
    pdf(file=paste(dir, "/", "normalization_comparison", ".pdf", sep=""), paper="a4r", width=11.69, height=8.27)
    my.normalization.comparison(kIn1.subset, normalizationMethods)
    dev.off()
    
    #hit detection using +- k MAD or SD
    pdf(file=paste(dir, "/", "significant_hit_count", ".pdf", sep=""), paper="a4r", width=11.69, height=8.27)
    qplot(Transfection, data=subset(significantHitList, !is.na(Transfection)), main="Pilot XV: Significant hits", stat="bin", fill=normalization) + scale_fill_brewer()
    dev.off()
  }
  
  return(significantHitList)
}

#calclate my own zscore
#use only 80 inner wells for computation since the outer columns are edge-biased and contain positive controls
#use negative controls on the outer wells for computation of mean and sd, median and mad
my.zscore <- function(plateColumn, CTblue, onlyInner=T, robust=F, na.rm=T)
{
  newdf <- data.frame(col=plateColumn, signal=CTblue)
  
  if(onlyInner)
    sub <- subset(newdf, col > 1 & col < max(col))$signal
  else
    sub <- newdf$signal
  
  if(!robust)
    return((CTblue - mean(sub, na.rm=na.rm))/ sd(sub,na.rm=na.rm))
  else
    return((CTblue - median(sub, na.rm=na.rm))/ mad(sub, na.rm=na.rm))
}


#box plots
my.boxplot <- function(kIn1, signalType, plotCV=T)
{
  title <- paste("Pilot XV", signalType ,": Boxplots (numbers correspond to CV)")
  p <- qplot(Transfection, signal, data=kIn1, ylab=signalType, main=title) + geom_boxplot(aes(fill = InnerOuter)) 
  p <- p + facet_wrap(~cellLine) + opts(axis.text.x=theme_text(angle=45))
  if(plotCV)
  {
    p <- p + stat_summary(fun.data=give.cv, geom="text")
  }
  print(p)
}

#time-course per transfection
my.plateComparison <- function(kIn1, signalType)
{
p <- qplot(plate, signal, data=kIn1, colour=TopBottom, main=paste("Pilot XV",signalType,": Plate Comparison per Transfection"), shape=LeftRight) + facet_wrap(~Transfection) + stat_smooth(aes(group=1))
print(p)
}

#outliers
my.outliers <- function(kIn1.subset, method, margin=2, onlyLower=F)
{
  #subset outliers, e.g. for labels in heatmap
  kIn1.mean <- mean(kIn1.subset$signal, na.rm=T)
  kIn1.median <- median(kIn1.subset$signal, na.rm=T)
  kIn1.mad <- mad(kIn1.subset$signal, na.rm=T)
  kIn1.sd <- sd(kIn1.subset$signal, na.rm=T)
  
  upper_margin <- margin
  
  if(onlyLower)
  {
    upper_margin <- 1000
  }
  
  if(method=="sd")
    kIn1.outliers <- subset(kIn1.subset, plateColumn > 1 & plateColumn < 12 & ((signal > kIn1.mean + upper_margin * kIn1.sd) | (signal < kIn1.mean - margin * kIn1.sd))) 
  else if(method=="mad")
    kIn1.outliers <- subset(kIn1.subset, plateColumn > 1 & plateColumn < 12 & ((signal > kIn1.median + upper_margin * kIn1.mad) | (signal < kIn1.median - margin * kIn1.mad)))
  else if(method=="none") 
    kIn1.outliers <- subset(kIn1.subset, plateColumn > 1 & plateColumn < 12 & (signal > upper_margin | signal < - margin))
  
  return(kIn1.outliers)
}

#heatmap
my.heatmap <- function(kIn1.subset, signalType, margin=2, ncol=3, method="sd")
{
  kIn1.outliers <- my.outliers(kIn1.subset, method, margin)
  
  require(ggplot2);
  title <- paste("Pilot XV", signalType ,": Heatmaps (samples with median +-",margin,"SD are labeled)")
  p2 <- qplot(plateColumn, plateRow, data=kIn1.subset, main=title, xlab="column", ylab="row")
  p2 <- p2 + geom_tile(line=0, aes(fill = signal));
  #p2 <- p2 + scale_fill_gradient2(midpoint=mean(kIn$log10, na.rm=T));
  p2 <- p2 + scale_fill_gradient(low = "yellow", high = "red", CONTOUR=FALSE, name=paste(signalType));
  p2 <- p2 + facet_wrap(~plate, ncol=ncol);
  p2 <- p2 + scale_x_continuous(expand=c(0,0));
  p2 <- p2 + geom_text(data=kIn1.outliers, aes(label=Transfection), vjust=-1.5, hjust=.7)
  p2 <- p2 + geom_segment(data=kIn1.outliers, aes(x=plateColumn-0.5, y=plateRow-0.5, xend=plateColumn, yend=plateRow), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
  p2 <- p2 + opts(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank(), panel.margin=unit(0.1, "lines"), panel.margin=unit(0, "lines"), plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"), 
  plot.title=theme_text(size=18), strip.background=theme_rect(fill="grey90", colour="grey50"))
  print(p2 + scale_y_reverse(expand=c(0,0)));
}

#well to well comparison on all plates
my.well.to.well <- function(kIn, signalType)
{
  require(ggplot2)
  f <- qplot(plate, signal, data=kIn, colour=Control, ylab=signalType, xlab="plate", main=paste("Pilot XV", signalType,": per well comparison")) + facet_grid(plateRow~plateColumn) + stat_smooth(aes(group=1))
  g <- f + geom_text(aes(5, 0, label=Transfection), vjust=-1, size=3) + opts(axis.text.x=theme_text(angle=-45, size=5))
  print(g)
}

#comparison of different normalization techniques
my.normalization.comparison <- function(kIn1, methods = c("CTblue", "Bscore", "zscore", "zscore_whole_plate", "poc", "npi"))
{
  normalized <- kIn1[,c("Transfection", "plate", methods)]
  normalized.melted <- melt(normalized, id.vars=c("plate", "Transfection"))
  normalized.melted$plate <- as.factor(normalized.melted$plate)
  p <- qplot(plate, value, data=normalized.melted) + facet_grid(variable~Transfection, scales="free") + stat_smooth(aes(group="1"))
  print(p)
}

#qqplot
my.qqplot <- function(kIn1, signalType, includeControls=F, margin=2)
{
  if(includeControls == F)
    data <- subset(kIn1, plateColumn > 1 & plateColumn < max(plateColumn))
  else
    data <- kIn1
    
  sorted <- arrange(data, signal)
  
  sorted$x <- 1:length(sorted$signal)/length(sorted$signal)
  significant <- subset(sorted, signal < -margin | signal > margin)
  
  p <- qplot(x, signal, data=sorted, xlim=c(0,max(sorted$signal)), ylim=c(0,max(sorted$signal))) + stat_abline(intercept=0,slope=1, col="red")
    p <- p + opts(title=paste("Pilot XV: QQ-Plot of", signalType))
    p <- p + scale_x_continuous(name="expected")
    p <- p + scale_y_continuous(name="observed")
    p <- p + geom_text(data=significant, aes(x=x, y=signal, label=Transfection, colour=factor(plate)))
    print(p)
}

#helper methods
give.cv <- function(x){
return(c(y = mean(x), label = (round(100 * sd(x) / mean(x)))))
}

give.sd <- function(x){
return(c(y = mean(x), label = round(sd(x))))
}

give.n <- function(x){
return(c(y = mean(x), label = (round(100 * sd(x) / mean(x)))))
}