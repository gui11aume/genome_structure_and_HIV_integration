print("This script calculates mean of how many times out did each real integration site score (for value of epigenomic feature in some bin) better than 10 random controls. Mean is then calculated over entire dataset.
IS_withaddedvalues_075.Robj contains dataset IS with value for epigenomic sites in all bins. IS_controls_075.Robj contins same values for all matched controls.
")

library(parallel)
library(stringr)
library(data.table)
library(pheatmap)
load("IS_withaddedvalues_075.Robj")
load("IS_controls_075.Robj")
compare <- function(x){
     ifelse(dt$RealValue>x,1,ifelse(dt$RealValue==x,1/2,0))
}


make_plot <- function(dat, group) {
     d <- as.data.frame(dat[,1:(ncol(dat)-1)])
     rownames(d) <- factor(dat$ChIP)
     d <- d[c("H3K27Ac", "H3K4me1", "BRD4", "Med1_RA", "H3K36me3", "H4K20me1","H3K4me3", "H3K27me3", "H3K9me2"),]
     print(
          pheatmap(d, main=group, cluster_rows=FALSE,cellwidth=40, cellheight=40)
     )
     
}    

IS_control <- IS_controls_075
IS <- IS_075

listNames <- colnames(IS_075)
epiNames <- rownames(IS_075)

GetThePhi <- function(epiNameWithBinsize, listName){
    require(data.table)
    dt <- as.data.table(IS_control[epiNameWithBinsize,grep(paste(listName,"_", sep=""), colnames(IS_control))])
    dt[,RealValue := as.data.table(IS[epiNameWithBinsize,listName])]
    a <- dt[,1:10]
    rowMeans(ifelse(a<dt$RealValue, 1,ifelse(a==dt$RealValue,1/2,0)))
}
###ROC by brady-adapted by Kuzman
ROC.MRC.pvalcalculation <- function(myphi_075)
{   
    names(myphi_075) <- colnames(IS_075)
    origin.levels <- colnames(IS_075)
    nvars <- ncol(myphi_075[[1]])
    phi.list <- myphi_075    
    rocz <- sapply(phi.list,colMeans)
    ## inflate the variance by 1e-10 to avoid diff/0.0 in Stats 
    roczVar <- lapply(phi.list,function(x) 1e-10*diag(ncol(x)) + var(x)/nrow(x))
    nullStats <- (rocz-0.5)^2/sapply(roczVar,diag)
    nullPvals <- pchisq(nullStats,df=1,lower.tail=FALSE)
    variableDiffs <-
        do.call(rbind,combn(nvars,2,function(x) rocz[x[1],]-rocz[x[2],],
                            simplify=FALSE))
    variableDVars <-
        sapply(roczVar,
               function(x) combn(nvars,2,
                                 function(y) sum(x[y,y]*c(1,-1,-1,1))))
    variableDStats <- variableDiffs^2/variableDVars
    variablePvals <- pchisq(variableDStats,df=1,lower.tail=FALSE)
    attr(variablePvals,"whichRow") <- combn(nvars,2)
    
    #        
    if (length(origin.levels)>1){
        originDVars <- combn(roczVar,2,function(x) diag(x[[1]])+diag(x[[2]]))
        originDiffs <- combn(origin.levels,2,function(x) rocz[,x[1]]-rocz[,x[2]])
        originStats <- originDiffs^2/originDVars
        originPvals <- pchisq(originStats,df=1,lower.tail=FALSE)
        attr(originPvals,"whichCol") <- combn(length(origin.levels),2)
        rownames(originPvals) <- rownames(rocz)
    } else {
        originPvals <- matrix(NA,nrow=nvars,ncol=1L)
    }
    
    DrawNP <- ifelse(nullPvals>0.05,"",ifelse(nullPvals>0.01,"*",ifelse(nullPvals>0.001,"**","***")))
    pvals <- list(vp=variablePvals,np=nullPvals, 
                  DrawNP=DrawNP, 
                  op=originPvals)
    list(ROC=rocz,var=roczVar,pvalues=pvals)
}

myphi_075 <- mclapply(listNames, function(listName){
    sapply(epiNames, function(epiName){
        GetThePhi(epiName, listName)
    })
}
#,mc.cores=14) # use if multiple cores enabled
,mc.cores=1)


roc<-ROC.MRC.pvalcalculation(myphi_075)

allValues_075 <- roc$ROC
colnames(allValues_075) <- colnames(IS_075)
drawMat <- roc$pvalues$DrawNP

dt <- as.data.table(allValues_075)
dt$SectionName <- rownames(allValues_075)
dt$ChIP <- substr(str_extract(dt$SectionName, "[A-Z].*_\\d"), 1,nchar(str_extract(dt$SectionName, "[A-Z].*_\\d"))-2)
dt$binSize<-substr(str_extract(dt$SectionName, "_\\d.*"), 2,nchar(str_extract(dt$SectionName, "_\\d.*")))

DRAWdt <- as.data.table(drawMat)
DRAWdt$SectionName <- rownames(allValues_075)
DRAWdt$ChIP <- substr(str_extract(dt$SectionName, "[A-Z].*_\\d"), 1,nchar(str_extract(dt$SectionName, "[A-Z].*_\\d"))-2)
DRAWdt$binSize<-substr(str_extract(dt$SectionName, "_\\d.*"), 2,nchar(str_extract(dt$SectionName, "_\\d.*")))

pdf("ROC_mainFigure_bybinsizes.pdf", width=7.75, height=9)
epigenNamesChosen <- c("BRD4", "H3K27Ac", "H3K27me3","H3K36me3","H3K9me2","Med1_RA", "H3K4me3", "H4K20me1", "H3K4me1")
dt[ChIP%in%epigenNamesChosen, make_plot(.SD, .BY), by=binSize, .SDcols=colnames(dt)[c(3:9,15,17)]]
dev.off()

epigenNamesChosen <- c("BRD4", "H3K27Ac", "H3K27me3","H3K36me3","H3K9me2","Med1_RA", "H3K4me3", "H4K20me1", "H3K4me1")
pdf("ROC_SupplementaryFigure_bybinsizes.pdf", width=3.75, height=9)
dt[ChIP%in%epigenNamesChosen, make_plot(.SD, .BY), by=binSize, .SDcols=colnames(dt)[c(12:14,17)]]
dev.off()

pdf("Figure1B.pdf", width=7.75, height=9)
make_plot(dt[(ChIP%in%epigenNamesChosen)& binSize==1000, c(3:9,15,17)], "1000")
dev.off()     

pdf("SupplementaryFigure1E.pdf", width=7.75, height=9)
dt[(ChIP%in%epigenNamesChosen)& binSize==1000, make_plot(.SD, .BY), by=binSize, .SDcols=colnames(dt)[c(12:14,17)]]
dev.off()

dff <- as.data.frame(dt[ChIP%in%epigenNamesChosen, c(3:9,12:15), with=FALSE])
rownames(dff) <- dt[ChIP%in%epigenNamesChosen,SectionName]
DRAWdff <- as.data.frame(DRAWdt[ChIP%in%epigenNamesChosen, c(3:9,12:15), with=FALSE])
rownames(DRAWdff) <- DRAWdt[ChIP%in%epigenNamesChosen]$SectionName
pdf("ROC_withasterisks.pdf", 8,18)
pheatmap(dff, cluster_rows=FALSE, display_numbers=DRAWdff, clustering_method="median")
dev.off()


pdf("ROC_allValues.pdf", 10,20)
pheatmap(allValues_075, cluster_rows=F, display_numbers=T)
dev.off()

allVars <- ls()[ls()!="roc"]

print(make_plot(dt[(ChIP%in%epigenNamesChosen)& binSize==1000, c(3:9,15,17)], "1000"))
print(dt[(ChIP%in%epigenNamesChosen)& binSize==1000, make_plot(.SD, .BY), by=binSize, .SDcols=colnames(dt)[c(12:14,17)]])
rm(list=allVars)
return(roc)
