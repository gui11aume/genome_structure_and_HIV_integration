require(ggplot2)
# plot for figure - "saturation"
# number of new genes found by adding another list, lists are sorted by 
# number of integration sites, descending.
library(GenomicRanges)
load("is.Robj")
load("genidf.Robj")
geni <- GRanges(genidf$seqnames, IRanges(genidf$start, genidf$end), Names=genidf$Names)
is <- is[c("Lucic", "Maldarelli", "Cohn", "Han", "Ikeda","Wagner", "Brady", "Kok")]
is <- is[order(sapply(is, length), decreasing = T)]
is

NumberOfUniqueGenesFound <- c()

for(i in 1:length(is)){
     namesOfdataset_i <- unique(as.data.frame(geni[countOverlaps(geni,is[[i]], ignore.strand=T)>0])$Names)
     NumberOfUniqueGenesFound[[i+1]] <- as.character(unique(c(NumberOfUniqueGenesFound[[i]],namesOfdataset_i)))
}

numberofuniquegenesHit <- sapply(NumberOfUniqueGenesFound, length)
numberOfISAddedToList <- cumsum(c(0,sapply(is,length)))
ddf <- data.frame(numberOfISAddedToList,numberofuniquegenesHit)
plotSaturation <- ggplot(ddf, aes(numberOfISAddedToList,numberofuniquegenesHit))+
     geom_smooth(method='lm', size=0.75)+
     geom_line(stat='identity', size=1.5, col="red")+
     theme_bw() +
     coord_cartesian(xlim=c(0,max(numberOfISAddedToList)), ylim=c(0,max(numberofuniquegenesHit)))+
     ylab("Number of genes hit by integration")+xlab("Number of integration sites")

lm_summary <- summary(lm(ddf$numberofuniquegenesHit~ddf$numberOfISAddedToList))

print(plotSaturation)

ggsave(plotSaturation, file="SupplementaryFigure1B.pdf", width = 7, height = 5)
