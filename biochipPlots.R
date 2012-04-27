p <- qplot(column, row, data=kIn)
p <- p + geom_tile(data= subset(kIn, !is.na(`F635 Mean - B635`)), aes(fill = `F635 Mean - B635`), colour = "black") 
p <- p + geom_tile(data= subset(kIn, is.na(`F635 Mean - B635`)), aes(colour = NA), fill="steelblue")
p <- p + coord_cartesian(ylim=c(14.5,0.5))
p <- p + scale_fill_gradient(low = "black", high = "red")
p <- p + facet_wrap(~block, ncol=11)
p <- p + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

duplicateRows <- function(kIn)
{
  for(rowBlocks in 0:3)
  {
    currentRow <- rowBlocks * 16 + 1
    kIn[(currentRow + 7):(currentRow + 13), 2:23] <- kIn[currentRow:(currentRow+6), 2:23]
  }
  return(kIn)
}

library(ggplot2);

kIn$log2 <- log2(kIn$normalized)
p4 <- qplot(column, row, data=kIn, main="log2(CD44 / GAPDH) Normalization")
p4 <- p4 + geom_tile(data= subset(kIn, !is.na(`normalized`)), line=0, aes(fill = `log2`), colour = "black") 
p4 <- p4 + geom_tile(data= subset(kIn, is.na(`normalized`)), line=0, aes(colour = NA),  fill="black")
p4 <- p4 + coord_cartesian(ylim=c(16.5,0.5));
p4 <- p4 + scale_fill_gradient2(midpoint=0, contour = TRUE);
p4 <- p4 + facet_wrap(~block, ncol=11);
p4 <- p4 + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0));
p4 <- p4 +opts(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank(), panel.margin=unit(0.1, "lines"), panel.margin=unit(0, "lines"), plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"), 
             plot.title=theme_text(size=18), strip.background=theme_rect(fill="grey90", 
colour="grey50"))
print(p4 + ylim(c(16.5,0.5)));

kIn.filtered <- kIn[!is.na(kIn$replicate),]
kIn.filtered <- kIn.filtered[!is.na(kIn.filtered[["F635 Mean - B635"]]),]
kIn.filtered$log2 <- log2(kIn.filtered[["F635 Mean - B635"]])
kIn.filtered$depositions <- as.factor(kIn.filtered$depositions)
p5 <- qplot(value, log2, data=kIn.filtered, main="Comparison FG-BG values")
p5 <- p5 + geom_boxplot(aes(fill = factor(depositions)))
p5 <- p5 + facet_wrap(~staining)
p5 <- p5 + opts(axis.text.x=theme_text(angle=-45))
p5 <- p5 + labs(fill = "Depositions")


theme_map <- function(size=12) 
{ 
    o = list(axis.line=theme_blank(), 
             axis.text.x=theme_blank(), 
             axis.text.y=theme_blank(), 
             axis.ticks=theme_blank(), 
             axis.ticks.length=unit(0.3, "lines"), 
             axis.ticks.margin=unit(0.5, "lines"), 
             axis.title.x=theme_blank(), 
             axis.title.y=theme_blank(), 
             legend.background=theme_rect(fill="white", colour=NA), 
             legend.key=theme_rect(colour="white"), 
             legend.key.size=unit(1.2, "lines"), 
             legend.position="right", 
             legend.text=theme_text(size=size*0.8), 
             legend.title=theme_text(size=size*0.8, face="bold", 
hjust=0), 
             panel.background=theme_blank(), 
             panel.border=theme_blank(), 
             panel.grid.major=theme_blank(), 
             panel.grid.minor=theme_blank(), 
             panel.margin=unit(0, "lines"), 
             plot.background=theme_blank(), 
             plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"), 
             plot.title=theme_text(size=size*1.2), 
             strip.background=theme_rect(fill="grey90", 
colour="grey50"), 
             strip.text.x=theme_text(size=size*0.8), 
             strip.text.y=theme_text(size=size*0.8, angle=-90)) 
    return(structure(o, class="options")) 
    
    
rows <- dcast(test, Column ~ Block, fun.aggregate=mean)[1:12,2:45]
rows.cor <- cor(rows)
corrgram(rows.cor, order=FALSE, main="Correlation of columns (row mean per block)", upper.panel=panel.pie)

cols <- dcast(kIn, Row ~ Block, fun.aggregate=mean)[1:14,2:45]
cols <- as.matrix(cols)
cols[is.nan(cols)] <- 0
cols.cor <- cor(cols)
corrgram(cols.cor, order=FALSE, main="Correlation of rows (column mean per block)", upper.panel=panel.pie)   
    
my <- as.data.frame(cbind(seq(1:44), apply(cols.cor, 1, mean)))
colnames(my) <- c("Block", "Correlation")
q<-qplot(data=my, x=Block, y=Correlation, main="Average Correlation of Blocks")
q<- q + geom_point(data=subset(my, Correlation < (1-2*sd(my$Correlation))), colour="red")
q + geom_text(aes(label=Block), data=subset(my, Correlation < (1-2*sd(my$Correlation))), vjust=1.2)
my<-cbind(my, RowShift=0)
my[subset(my, Correlation < (1-2*sd(my$Correlation)))$Block,] <- 1
    
myValues <- kIn[!is.na(kIn[["signal"]]),c("label", "signal", "cellEquivalent", "depositions")]
myValues$cellEquivalent[is.na(myValues$cellEquivalent)] <- "_1"
myValues$cellEquivalent <- factor(myValues$cellEquivalent, levels=c("_1", "_2", "_4", "_8", "_16"), ordered=T)
myValues<-cast(myValues, label + depositions ~cellEquivalent, value="signal", add.missing=TRUE, fun.aggregate="median", na.rm=T)
data.dilutes <- as.matrix(myValues[,3:ncol(myValues)] )  
    
serialDilutionParams <- plot.dilution.series(D0=2, data.dilutes=data.dilutes, sensible.max=1.e4)

D0 <- 2 # preset dilution factor, known from experiments
    
a=serialDilutionParams$a; c=serialDilutionParams$c; D=serialDilutionParams$D; 
d.a=serialDilutionParams$d.a; d.c=serialDilutionParams$d.c; d.D=serialDilutionParams$d.D; 
    
serialDilutionResult <- protein.con(D0=D0, D=D,c=c,a=a,d.a=d.a, d.D=d.D, d.c=d.c,data.dilutes=data.dilutes)

X.result = serialDilutionResult
x.saturation.level=   D0^(3)/((1/( M/r - a)- 1/(M-a)))^(log(D0)/log(D)) 
x.nodetection.level = D0^(0)/((1/( r*a - a)- 1/(M-a)))^(log(D0)/log(D))

#plot(log2(X.result[,1]),log2(X.result[,2]))
plot( log2(X.result[,1]) , X.result[,2]/X.result[,1],xlab='log2(Estimated Con)', ylab='CV',ylim=c(0,0.2),xlim=c(4,20))
abline(v=log2(x.saturation.level),lty=4,col='green')
abline(v=log2(x.nodetection.level),lty=4,col='blue')
    
#making plot Signal vs estimated concentration
#X.result[,1] =x.t
plot(log2(X.result[,1]), log2(myValues[,3]),col='red',pch='+', xlab='log2(Estimated Con)', ylab='log2(Signal)')
points(log2(X.result[,1]/2), log2(myValues[,4]),col='orange',pch='+')
points(log2(X.result[,1]/4), log2(myValues[,5]),col='green',pch='+')
points(log2(X.result[,1]/8), log2(myValues[,6]),col='blue',pch='+')
  
    