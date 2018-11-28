library(ggplot2)
library(tidyr)
library(GenomicRanges)
load("is.Robj")
nLists <- 9
load("genidf.Robj")
geni <- GRanges(genidf$seqnames, IRanges(genidf$start, genidf$end), strand=genidf$strand,
                Names=genidf$Names)
# Finding genes that overlap with unique IS 
GeneNames <- sapply(is, function(x){
     unique(geni[countOverlaps(geni, x, ignore.strand=T)>0,]$Names)
})
# I will add two more sets of random genes pooled from union of all hit CD4HIV genes
# - one with length as the smallest IS data set
# and the other as the largest 
set.seed(2163)
smallestRandom <- replicate(1000,sample(unique(geni$Names),length(GeneNames$Han)-1))
largestRandom <- replicate(1000,sample(unique(geni$Names),length(GeneNames$Cohn)+1))

#This is just a dummy set, values are replaced later with mean + sd of 1000 repetitions
GeneNames$smallestRandom <- smallestRandom[,1]
GeneNames$largestRandom <- largestRandom[,1]


sortingorderbysize <- order(sapply(GeneNames,length))
# what percentage of genes found in a list is shared with other lists
PercentageOfOverlapWithLists <- sapply(1:nLists, function(i){
     mean(GeneNames[[i]] %in% unique(unlist(GeneNames[!(1:nLists)%in%i])))
})

#calculation Of percent overlaps for smallest and largest random genes (X1000)
PercentageOfOverlapWithLists_smallest <- sapply(1:ncol(smallestRandom), function(i){
     mean(smallestRandom[,i] %in% unique(unlist(GeneNames[!(1:nLists)%in%i])))
})
PercentageOfOverlapWithLists_largest <- sapply(1:ncol(largestRandom), function(i){
     mean(largestRandom[,i] %in% unique(unlist(GeneNames[!(1:nLists)%in%i])))
})
PercentageOfOverlapWithLists[nLists+1] <- mean(PercentageOfOverlapWithLists_smallest)
PercentageOfOverlapWithLists[nLists+2] <- mean(PercentageOfOverlapWithLists_largest)
# what is the percentage of shared genes between pairs of lists
PointPercentageOfOverlapBetweenLists <- sapply(1:length(GeneNames), function(i){
          sapply(1:length(GeneNames), function(j){
          mean(GeneNames[[i]] %in% GeneNames[[j]] )
          })
})
colnames(PointPercentageOfOverlapBetweenLists) <- rownames(PointPercentageOfOverlapBetweenLists) <-names(GeneNames)

# # version 1: mini pies 
PPercOfOvl_forPlots <- PointPercentageOfOverlapBetweenLists[sortingorderbysize,sortingorderbysize]
# SupplementaryFigure1_version1 <- corrplot(PPercOfOvl_forPlots, method = "pie")

# version 2: boxplot of pairwise overlaps
PPercOfOvl_forPlots <- as.data.frame(PPercOfOvl_forPlots)
PPercOfOvl_forPlots$OriginalList <- rownames(PPercOfOvl_forPlots) 
ppercovl_long <- gather(PPercOfOvl_forPlots, List1, PercentageOverlap, smallestRandom:largestRandom)
ppercovl_long <- ppercovl_long[-which( ppercovl_long$OriginalList==ppercovl_long$List1     ),]
ppercovl_long$OriginalList <- factor(ppercovl_long$OriginalList, levels=rownames(PPercOfOvl_forPlots))
ppercovl_long$List1 <- factor(ppercovl_long$List1, levels=rownames(PPercOfOvl_forPlots))

SupplementaryFigure1_version2 <- ggplot(ppercovl_long, aes(List1, PercentageOverlap))+
     ylim(0,1)+
     theme_bw()+
     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     ylab("       Proportion of genes that 
          overlap with other lists")+
     xlab("")+
     stat_boxplot(geom="errorbar")+
     geom_boxplot()

# version 3: barplot of overlaps with all other lists
PercentageOfOverlap <- data.frame(Names=names(GeneNames), 
                                  Overlap = PercentageOfOverlapWithLists)

PercentageOfOverlap$Names <- factor(as.character(PercentageOfOverlap$Names),
                                    levels=names(GeneNames)[sortingorderbysize])
PercentageOfOverlap$ymin <- PercentageOfOverlap$Overlap
PercentageOfOverlap$ymax <- PercentageOfOverlap$Overlap
PercentageOfOverlap$ymax[nLists+2] <- quantile(PercentageOfOverlapWithLists_largest,0.95)
PercentageOfOverlap$ymax[nLists+1] <- quantile(PercentageOfOverlapWithLists_smallest,0.95)

SupplementaryFigure1_version3 <- ggplot(PercentageOfOverlap, aes(Names,Overlap))+
     geom_bar(stat="identity", fill="white",color="black")+
     ylim(0,1)+
     theme_light()+
     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     ylab("Proportion of genes that overlap with 
     any of the other lists")+
     xlab("")+
     geom_errorbar(aes(ymin = ymin,ymax = ymax))
message("Code generates 3 plots:
SupplementaryFigure1_version1: pairwise percentage of shared genes
SupplementaryFigure1_version2: boxplot of proportion of genes shared with each of other list 
SupplementaryFigure1_version3: barplot of proportionof genes in a list that is shared with at 
least one other list 
        Only version 3 is used in the paper as suplementary figure3:")

print(SupplementaryFigure1_version3)
ggsave(SupplementaryFigure1_version3, file="SupplementaryFigure3.pdf", width=7, height=5)
