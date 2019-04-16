load("is.Robj")
load("GRCh37AllGenes.Robj")
library(xlsx)
is$InVitro <- c(is$Lucic, is$Brady)
is$Patients <- c(is$Maldarelli, is$Cohn, is$Han, is$Ikeda, is$Kok, is$Wagner)    
is$`Total CD4 HIV` <- unique(c(is$Maldarelli, is$Cohn, is$Han, is$Ikeda, is$Kok, is$Wagner, is$Lucic, is$Brady))

TableM1<- sapply(is, function(x)
{
     data.frame(
          totalNumberOfIS = length(x),
          numberOfUniqueSitesInGenes = sum(countOverlaps(x,geni, ignore.strand=T)>0),
          percentOfUniqueSitesInGenes = round(mean(countOverlaps(x,geni, ignore.strand=T)>0),2),
          numberOfGenes = length(geni),
          numberOfGenesWithIS = sum(countOverlaps(geni,x, ignore.strand=T)>0),
          percentOfGenesWithIS = round(mean(countOverlaps(geni,x, ignore.strand=T)>0),2),
          numberOfGenesWithoutIS = sum(countOverlaps(geni,x, ignore.strand=T)==0)
     )
})
rownames(TableM1) <- c("Total number of unique IS",
                             "# of unique IS in genes",
                             "% of unique IS in genes",
                             "Total number of genes",
                             "# of genes with IS",
                             "% of genes with IS",
                             "# of genes without IS")

message("TableM1 contains all basic gene statistics")
write.xlsx(TableM1, file="Additional file 1.xlsx")


