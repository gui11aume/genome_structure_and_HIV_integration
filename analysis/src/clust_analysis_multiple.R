nsuppressMessages(library(GenomicRanges))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(Sushi))

## Params
hic.binsize  = 5e3
expr.cutoff  = 0.4
enh.radius   = 2.5e3
prom.radius  = 2.5e3
hp.binsize   = 1e5
hphr.binsize = 1e4

## Files
fig.dir      = 'figures/'
cmarks.files = 'chip_datasets/reduced/*.bed'
cmarks.type  = 'reduced' # all or reduced
cluster.list = list('oe_ice'='../clusters/global_clusters.txt')
global.clust.names = list('oe_ice'=list('0'='A2','1'='B1','2'='B3','3'='A1','4'='B2'),
                          'oe_none'=list('0'='B3','1'='A1','2'='B1','3'='B2','4'='A2'),
                          'raw3_ice'=list('0'='A2','1'='B3','2'='A1','3'='B2','4'='B1'),
                          'raw3_none'=list('0'='A1','1'='B3','2'='B2','3'='B1','4'='A2'))
hiv.datasets = 'hiv_datasets/*.integ'
hiv.expr     = 'expr_datasets/bhive_1_expr.txt'
gene.expr    = 'expr_datasets/jurkat_rnaseq.tpm'
chrsz.data   = 'annotations/hg19_chromsizes.txt'
enh.data     = 'annotations/H3K27ac.01'
hp.ctrl.hk36 = 'annotations/chip_peaks/H3K36me3.broadPeak'
hp.ctrl.pol2 = 'annotations/chip_peaks/PolII.narrowPeak'
hp.samp.size = 200
hp.samp.file = 'annotations/hp_ctrl.txt'
hiv.gen.file = 'annotations/hiv_gen.txt'
sup.enh.file = 'annotations/Jurkat_SE_all.bed'
hp.sen.file  = 'annotations/hp_sen.txt'
hp.file      = 'annotations/hiv_hotspots.txt'
sen.pcs.file = 'annotations/super_enh_ranges'
cool.file    = 'cool/WT_hg19_5k_q10.cool'
hiv.gen.len  = 1e3
hiv.gen.samp = 200
rig.gen.tss  = 1e3 ## Offset (upstream and downstream) from TSS
rig.gen.tes  = 1e3 ## Offset (upstream and downstream) from TES
rig.gen.step = 50  ## Bin step
rig.gen.bin  = 100 ## Bin width

## Analysis setup
pattern_sizes_pie    = TRUE
pattern_category_pie = TRUE
intrachromosomal_chromatin_heatmaps = TRUE
interchromosomal_chromatin_heatmap  = TRUE
interchromosomal_hiv_heatmap        = TRUE
interchromosomal_chromatin_heatmaps_light = TRUE
interchromosomal_contact_tests = TRUE


## Function handlers
setdiff_nw      <- function(...) { suppressWarnings(setdiff(...)) };
findOverlaps_nw <- function(...) { suppressWarnings(findOverlaps(...)) };

## Helper functions
MaxTable <- function(InVec) {
    if (!is.factor(InVec)) InVec <- factor(InVec)
    if (length(InVec) == 1)
        InVec
    t       <- table(InVec)
    maxVals <- names(t)[t == max(t)]
    if (length(maxVals) > 1)
        NA
    else
        maxVals
}


### File index
##
## 1. Read data.
## 1.1. Read chromatin marks.
## 1.2. Read HIV datasets.
## 1.3. Read annotations.
## 1.4. Read 3D patterns.
## 1.5. Read 3D pattern spectral info
##
## 2. Data structures.
## 2.1. Compute genomic categories.
## 2.2. 3D patterns.
## 2.3. Category background of 3D patterns
##
## 3. Signal distribution in 3D patterns.
## 3.1. Chromatin marks in 3D patterns.
## 3.2. HIV in 3D patterns.
##
## 4. HIV integration pattern
## 4.1. Integration density in Categories.
## 4.2. Integration density in 3D patterns.
## 4.3. Local integration density in categories (grouped and normalized by 3D pattern)
## 4.4. Local integration density in 3D patterns (grouped and normalized by Category)
## 4.5. HIV integration density (grouped by 3D pattern or category)
## 
## 5. Gene expression in 3D patterns
## 5.1. Endogenous expression in 3D patterns
## 5.2. HIV expression in 3D patterns
##
## 6. HIV hotspots
## 6.1. Hotspot occupancy
## 6.2. HIV expression in hotspots.
## 6.3. Hotspots and SuperEnhancers
## 6.4. Chromatin composition in Hotspots
## 6.5. Hotspots controls (PolII, H3K36me3)
##
## 8. RIG genes in Jurkat
## 8.1. Clustering of HIV-dense active genes
## 8.2. 3D pattern of integrations in HIV-dense gene
##
## 9. SuperEnhancers and HD genes.
## 9.1. Compute PCS between HD genes and SuperEnhancers (DBsuper)


##
## 1. Read data.
##
cat('[1. Read data]\n')


##
## 1.1. Read chromatin marks.
##
cat('[1.1. Read chromatin marks]\n')
cm.gr    = GRanges()
cm.names = c()

dataset.id = 1
for (f in Sys.glob(cmarks.files)) {
    ## Extract CM name from filename (text before '_')
    path     <- strsplit(f,'/')[[1]]
    fname    <- strsplit(path[length(path)],'.',fixed=T)[[1]][1]
    cm.name  <- paste(dataset.id,strsplit(fname,'_',fixed=T)[[1]][1],sep='.')
    cm.names <- c(cm.names,cm.name)
    ## Read file (only chr[1-9XY][0-9]*)
    d       <- subset(read.table(f), grepl("^chr[0-9XY][0-9]*$",V1))
    d$V1    <- factor(d$V1, unique(d$V1))
    cm.gr   <- c(cm.gr, GRanges(seqnames=Rle(d$V1), ranges=IRanges(start=d$V2,end=d$V3),
                                dataset= rep(dataset.id,nrow(d))))
    ## Verbose dataset info
    cat('[1.1] chromatin mark #',dataset.id,': ',cm.name,' (',nrow(d),' lines)\n',sep='')
    ## Increase loop counter
    dataset.id <- dataset.id + 1
}

##
## 1.2. Read HIV datasets.
##
cat('[1.3. Read HIV datasets]\n')

## 1.2.1. Insertion sites.
cat('[1.3.1. HIV insertion sites]\n')
hiv.gr    = GRanges()
hiv.names = c()

dataset.id = 1
for (f in Sys.glob(hiv.datasets)) {
    ## Extract dataset name from filename (text before '_')
    path      <- strsplit(f,'/')[[1]]
    fname     <- strsplit(path[length(path)],'.',fixed=T)[[1]][1]
    hiv.name  <- strsplit(fname,'_',fixed=T)[[1]][1]
    hiv.names <- c(hiv.names,hiv.name)
    ## Read file (only chr[1-9XY][0-9]*)
    d         <- subset(read.table(f), grepl("^chr[0-9XY][0-9]*$",V1))
    d$V1      <- factor(d$V1, unique(d$V1))
    data.gr   <- GRanges(seqnames = Rle(d$V1), ranges = IRanges(start=d$V2, width=1), strand = Rle(d$V3),
                         dataset      = rep(dataset.id, nrow(d)),
                         dataset.name = rep(hiv.name, nrow(d))
                         )
    hiv.gr    <- c(hiv.gr,data.gr)
    ## Verbose dataset info
    cat('[1.3.1] HIV dataset #',dataset.id,': ',hiv.name,' (',nrow(d),' insertion sites)\n',sep='')
    ## Increase loop counter
    dataset.id <- dataset.id + 1
}

## 1.2.2. HIV expression (BHIVE)
cat('[1.3.2. HIV expression data (BHIVE)]\n')
bhive    <- read.table(hiv.expr, header=T)
bhive.gr <- GRanges(Rle(bhive$chr), IRanges(start=bhive$locus, width=1),
                    strand=Rle(bhive$strand), brcd=bhive$brcd,
                    nread=bhive$reads, mapq=bhive$mapq, rep=bhive$rep, expr=bhive$exprscore)
## Verbose
cat('[1.3.2] hiv integrations with expression info: ',sum(!is.na(bhive$exprscore)),'\n')

##
## 1.3. Read annotations.
##
cat('[1.4. Read annotations]\n')

## 1.3.1. Reference genome.
cat('[1.4.1. Reference genome]\n')
genome         <- read.table(chrsz.data, as.is=TRUE, col.names=c('chr','size'))
genome.gr      <- GRanges(Rle(genome$chr), IRanges(start=1, width=genome$size))
genome$bins    <- ceiling(genome$size/hp.binsize)
genome$bins.hr <- ceiling(genome$size/hphr.binsize)
genome.bins    <- data.frame(chrom=rep(genome$chr,genome$bins), bin=unlist(lapply(genome$bins, seq))-1, count=0)
genome.bins.hr <- data.frame(chrom=rep(genome$chr,genome$bins.hr), bin.hr=unlist(lapply(genome$bins.hr, seq))-1, count=0)


## Verbose genome
cat('[1.4.1] chromosomes: ',length(genome.gr),'\n',sep='')
cat('[1.4.1] coverage: ',sum(as.numeric(sum(coverage(genome.gr)))),' nt\n',sep='')
cat('[1.4.1] hotspot bin size: ',hp.binsize,'\n',sep='')
cat('[1.4.1] hotspot bin count: ',sum(genome$bins),'\n',sep='')

## 1.3.2. Read expression/gene data.
cat('[1.4.2. Gene expression info]\n')
genes     <- subset(read.table(gene.expr,header=T), type=='protein_coding')
genes$chr <- paste('chr',genes$chr,sep='')
cutoff    <- quantile(genes$tpm, expr.cutoff)
genact    <- subset(genes, tpm > cutoff)
gensil    <- subset(genes, tpm <= cutoff)
## Verbose RNA-seq info.
cat('[1.4.2] genes: ',nrow(genes),'\n',sep='')
cat('[1.4.2] expression cutoff: ',cutoff,' tpm\n',sep='')
cat('[1.4.2] active genes: ',nrow(genact),'\n',sep='')
cat('[1.4.2] silent genes: ',nrow(gensil),'\n',sep='')

## 1.3.3 Read enhancer data.
cat('[1.4.3. Enhancer data]\n')
enh           <- read.table(enh.data, comm="#")
enh.gr        <- reduce(GRanges(Rle(enh$V1), IRanges(start=enh$V2, enh$V3)))
start(enh.gr) <- start(enh.gr) - enh.radius
end(enh.gr)   <- end(enh.gr) + enh.radius
cat('[1.4.3] enhancer bed file: ',enh.data,'\n',sep='')
cat('[1.4.3] enriched intervals: ',length(enh.gr),'\n',sep='')
cat('[1.4.3] coverage: ',sum(as.numeric(sum(coverage(enh.gr)))),' nt\n',sep='')

## 1.3.4 Read SuperEnhancer DB.
cat('[1.4.4. SuperEnhancer data]\n')
sen.data     <- read.table(sup.enh.file, header=F, col.names=c('chr','beg','end'))
sen.gr       <- reduce(GRanges(seqnames = sen.data$chr, ranges = IRanges(start=sen.data$beg, end=sen.data$end)))
sen.data     <- data.frame(chr = seqnames(sen.gr),
                           beg = start(sen.gr),
                           end = end(sen.gr),
                           id  = seq(1,length(sen.gr)),
                           cat = "SuperEnhancer",
                           val = 1)
cat('[1.4.4] enhancer bed file: ',sup.enh.file,'\n',sep='')
cat('[1.4.4] superenhancer intervals: ',length(sen.gr),'\n',sep='')
cat('[1.4.4] coverage: ',sum(as.numeric(sum(coverage(sen.gr)))),' nt\n',sep='')

##
## 1.4. Read 3D patterns
##
cat('[1.2. Read 3D patterns]\n')
clust.grs = list()

for (i in seq(1,length(cluster.list))) {
    clust.gr = GRanges()
    clust.name = names(cluster.list)[i]
    cat('[1.2] ', 'cluster: ', clust.name, '\n', sep='')
    cat('[1.2] ', 'file: ', cluster.list[[clust.name]], '\n', sep='')
    for (line in readLines(file(cluster.list[[clust.name]]))) {
        words      <- strsplit(line,split="\t")[[1]]
        ## Parse columns.
        chr        <- words[1]
        local.clu  <- as.integer(words[2])
        global.clu <- as.integer(words[3])
        clu.idx    <- as.integer(unlist(strsplit(words[4], split=",")))
        ## Assemble cluster GR object.
        clu.gr     <- GRanges(seqnames = chr, IRanges(start=clu.idx*hic.binsize, end=(clu.idx+1)*hic.binsize-1),
                              local.cluster  = rep(local.clu, length(clu.idx)),
                              global.cluster = rep(global.clust.names[[clust.name]][[words[3]]], length(clu.idx)),
                              chrom          = rep(chr, length(clu.idx))
                              )
        clust.gr   <- c(clust.gr,clu.gr)
    }
    ## Verbose clusters
    cat('[1.2] ', length(unique(clust.gr$chrom)), ' chromosomes\n', sep='')
    cat('[1.2] ', length(unique(clust.gr$local.cluster)), ' intrachromosomal 3D patterns\n', sep='')
    cat('[1.2] ', length(unique(clust.gr$global.cluster)), ' global 3D patterns\n', sep='')

    clust.grs[[clust.name]] <- clust.gr
}


##
## 2. Data structures.
##
cat('[2. Data structures]\n')


##
## 2.1. Compute genomic categories.
##
cat('[2.1. Genomic categories]\n')
genes.gr  <- GRanges(Rle(genes$chr), IRanges(start=genes$beg, end=genes$end),
                     strand=Rle(genes$strand), name=genes$gene_name, expr=genes$tpm, active=genes$tpm > cutoff)
actgen.gr <- GRanges(Rle(genact$chr), IRanges(start=genact$beg, end=genact$end),
                     strand=Rle(genact$strand), name=genact$gene_name, expr=genact$tpm)
silgen.gr <- GRanges(Rle(gensil$chr), IRanges(start=gensil$beg, end=gensil$end),
                     strand=Rle(gensil$strand), name=gensil$gene_name, expr=gensil$tpm)
actpro.gr <- promoters(actgen.gr, up=prom.radius, down=prom.radius)
acttss.gr <- resize(actgen.gr,1)
siltss.gr <- resize(silgen.gr,1)

## TODO: Why not new categories, like: gene+enhancer, promoter+enhacer, enhancer alone?
## TODO: Also for the promoter, would it make sense to split in internal, external? Or delete the category?

## Active genes prevail over silent ones.
silgen_na     <- setdiff_nw(silgen.gr, actgen.gr, ignore.strand = TRUE)
## Subtract promoter regions from genes.
actgen_np     <- setdiff_nw(actgen.gr, actpro.gr, ignore.strand = TRUE)
silgen_nanp   <- setdiff_nw(silgen_na, actpro.gr, ignore.strand = TRUE)
## Subtract enhancers from genes and promoters
actpro_ne     <- setdiff_nw(actpro.gr, enh.gr, ignore.strand = TRUE)
actgen_npne   <- setdiff_nw(actgen_np, enh.gr, ignore.strand = TRUE)
silgen_nanpne <- setdiff_nw(silgen_nanp, enh.gr, ignore.strand = TRUE)

## Create category ranges list
categ = list()
categ[['silgen']] = reduce(silgen_nanpne)
categ[['actgen']] = reduce(actgen_npne)
categ[['actpro']] = reduce(actpro_ne)
categ[['enh']]    = reduce(enh.gr)
categ[['int']]    = reduce(setdiff_nw(genome.gr,union(silgen_nanpne, union(actgen_npne, union(actpro_ne,enh.gr)) )))

## Verbose categories.
cat('[2.1] active genes: ',sum(as.numeric(sum(coverage(categ[['actgen']])))),' nt\n',sep='')
cat('[2.1] silent genes: ',sum(as.numeric(sum(coverage(categ[['silgen']])))),' nt\n',sep='')
cat('[2.1] active promoters:: ',sum(as.numeric(sum(coverage(categ[['actpro']])))),' nt\n',sep='')
cat('[2.1] enhancers: ',sum(as.numeric(sum(coverage(categ[['enh']])))),' nt\n',sep='')
cat('[2.1] intergenic: ',sum(as.numeric(sum(coverage(categ[['int']])))),' nt\n',sep='')

## Define category colors
categ.colors        <- c("#92C46DA0","#548B54", "#19485780", "#000000A0", "#C0955FE0")
names(categ.colors) <- c('actgen','actpro','enh','silgen','int')


##
## 2.2. 3D patterns.
##
cat('[2.2. 3D patterns]\n')

##
## 2.2.1. 3D pattern PCA, eig-scores and abscore.
##
cat('[2.2.1. Read 3D pattern spectral info]\n')

require(lattice)
eig.grs = list()
for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.file = cluster.list[[clust.name]]
    cat('[2.2.1] cluster name: ',clust.name,'\n',sep='')    
    cluster.eigs = paste(dirname(clust.file),'/*.weig',sep='')
    eig.gr = GRanges()
    for (f in Sys.glob(cluster.eigs)) {
        ## Read file (only chr[1-9XY][0-9]*)
        cat('[2.2.1] eig file: ',f,'\n',sep='')    
        d       <- read.table(f)
        if (f == paste(clust.file, '.weig', sep='')) {
            d$V1    <- as.factor(d$V1)
            df      <- data.frame(clu=d$V1, pca='PCA2', x=d$V2, y=d$V3)
            df      <- rbind(df, data.frame(clu=d$V1, pca='PCA3', x=d$V2, y=d$V4))
        } else {
            d$V2    <- as.factor(d$V2)
            df      <- data.frame(clu=d$V2, pca='PCA2', x=d$V3, y=d$V4)
            if (ncol(d) > 4) {
                df      <- rbind(df, data.frame(clu=d$V2, pca='PCA3', x=d$V3, y=d$V5))
            } else {
                df      <- rbind(df, data.frame(clu=d$V2, pca='PCA3', x=d$V3, y=d$V4))
            }
            ## Get chromosome from filename (clusters_chr##.weig)
            d$chr = strsplit(strsplit(basename(f),'\\.')[[1]][1],'_')[[1]][2]
            d$beg = d$V1*hic.binsize+1
            d$end = (d$V1+1)*hic.binsize
            d$eig = d$V3
            ## Correct sign with gene expression.
            tmp.gr = GRanges(seqnames = d$chr,
                             IRanges(start=d$beg, end=d$end),
                             eig = d$eig
                             )
            tmp.ov = findOverlaps(tmp.gr, genes.gr)
            if (cor(tmp.gr[tmp.ov@from]$eig, genes.gr[tmp.ov@to]$expr) < 0) {
                tmp.gr$eig = -tmp.gr$eig
            }
            ## append to eig list
            eig.gr = c(eig.gr,tmp.gr)
        }
        ## Plot scatterplot
        figure.fn = paste(f,'_',clust.name,'.pdf',sep='')
        cat('[2.2.1] Weighted eigenvector scatterplot -> ',figure.fn,'\n',sep='')
        pdf(figure.fn,useDingbats=F,width=6,height=10)
        print(ggplot(data=df, aes(x=x, y=y, color=clu)) +
              geom_point() +
              facet_wrap(~pca, scales = 'free_y', nrow=2, strip.position = 'left',
                         labeller = as_labeller(c(PCA2 = "PCA2", PCA3 = "PCA3"))) + 
              ylab(NULL) +
              xlab('PCA1')
              )
        dev.off()
    }
    eig.grs[[clust.name]] = eig.gr
}

abscore.grs = list()
for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.file = cluster.list[[clust.name]]
    cat('[2.2.1] cluster name: ',clust.name,'\n',sep='')    
    cluster.abs = paste(dirname(clust.file),'/*.abscore',sep='')
    abscore.gr = GRanges()
    for (f in Sys.glob(cluster.abs)) {
        ## Read file (only chr[1-9XY][0-9]*)
        cat('[2.2.1] A/B score file: ',f,'\n',sep='')
        d <- read.table(f)
        # Make genomic ranges with abscores.
        chr = strsplit(strsplit(basename(f),'\\.')[[1]][1],'_')[[1]][2]
        tmp.gr = GRanges(seqnames = chr,
                         IRanges(start = d$V1*hic.binsize+1,
                                 end = (d$V1+1)*hic.binsize),
                         abscore = d$V2)
        ## Correct sign with gene expression.
        tmp.ov = findOverlaps(tmp.gr, genes.gr)
        if (length(tmp.ov) > 0 && cor(tmp.gr[tmp.ov@from]$abscore, genes.gr[tmp.ov@to]$expr) < 0) {
            tmp.gr$abscore = -tmp.gr$abscore
        }
        ## append to eig list
        abscore.gr = c(abscore.gr,tmp.gr)
    }
    abscore.grs[[clust.name]] = abscore.gr
}

##
## 2.2.2. 3D pattern structures.
##


for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]
    ## Cluster count and size.
    clust.bins <- table(clust.gr$global.cluster)
    clust.sz   <- clust.bins * hic.binsize
    clust.cnt  <- length(clust.sz)

    ## Define 3D pattern colors
    clust.colors               <- rev(brewer.pal(clust.cnt,"Blues"))
                                        #names(clust.colors)        <- c("Act+","Act-","H3K27me","H3K9me","Lamin")
    names(clust.colors)        <- c("A1","A2","B1","B2","B3")

    ## 3D pattern sizes pie chart.
    if (pattern_sizes_pie) {
        figure.fn <- paste(fig.dir,'3dpattern_sizes_',clust.name,'.pdf',sep='')
        cat('[2.2] cluster sizes (',clust.name,') -> ',figure.fn,'\n',sep='')

        ## Pie chart, cluster sizes.
        pdf(figure.fn,useDingbats=F)
        labels = paste('Cluster ',names(clust.sz),'\n',clust.sz/1e6,' Mbp')
        pie(clust.bins,labels=labels,col=clust.colors[names(clust.bins)],main='3D pattern sizes')
        dev.off()
    }

}
##
## 2.3. Category background of 3D patterns
##
cat('[2.3 Category background of 3D patterns]\n')

clust.categs    <- list()
clust.categs.gr <- list()

## Compute 3D pattern / category intersect
for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]

    clust.categ    <- data.frame(clust=c(),cat=c(),span=c())
    clust.categ.gr <- GRanges()
    for (i in unique(clust.gr$global.cluster)) {
        clust.df <- data.frame(clust=c(),cat=c(),span=c())
        ## iterate over pattern categories
        for (idx in names(categ)) {
            gc             <- clust.gr[clust.gr$global.cluster == i]
            ## Intersect and annotate cluster and background category.
            clu.cat        <- intersect(gc,categ[[idx]],ignore.strand=T)
            clu.cat$cat    <- idx
            clu.cat$clust  <- i
            clust.categ.gr <- c(clust.categ.gr, clu.cat)
            ## Store intersect span.
            categ.df       <- data.frame(clust=toString(i),cat=idx,span=sum(as.numeric(sum(coverage(clu.cat)))))
            categ.df$cat   <- levels(droplevels(categ.df$cat))
            clust.df       <- rbind(clust.df,categ.df)
        }
        clust.categ <- rbind(clust.categ,clust.df)
    }

    ## Save to lists
    clust.categs[[clust.name]] <- clust.categ
    clust.categs.gr[[clust.name]] <- clust.categ.gr

    ## Plot background categories in pie charts.
    if (pattern_category_pie) {
        figure.fn = paste(fig.dir,'3dpattern_background_category_',clust.name,'.pdf',sep='')
        cat('[2.3] 3D pattern background category pie charts (',clust.name,') -> ',figure.fn,'\n',sep='')
        pdf(figure.fn,useDingbats=F,width=20,height=15)
        par(mfrow=c(2,ceiling(clust.cnt/2)))
        for (i in unique(clust.gr$global.cluster)) {
            clust.df <- clust.categ[clust.categ$clust == i,]
            pie(clust.df$span, labels=paste(clust.df$cat,'\n',round(clust.df$span/1e4)/1e2,' Mbp',sep=''),main=paste('cluster',i),col=categ.colors[clust.df$cat])
        }
        dev.off()
    }
}

##
## 2.4. Overlap between 3D pattern and A/B scores.
##

for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]
    eig.gr = eig.grs[[clust.name]]
    
    cat('[2.4] First eigenvector (',clust.name,')\n',sep='')
    eig.clust.ov <- findOverlaps(eig.gr, clust.gr, ignore.strand=TRUE)
    eig.gr$clust <- NA
    eig.gr$clust[eig.clust.ov@from] <- clust.gr[eig.clust.ov@to]$global.cluster
    eig.clu.df = data.frame(chr=seqnames(eig.gr), eig=eig.gr$eig, clust=eig.gr$clust)
    eig.clu.df = eig.clu.df[!is.na(eig.clu.df$eig) & !is.na(eig.clu.df$clust),]

    cat('[2.4] Distribution of first eigenvector by 3D pattern (',clust.name,')-> ',figure.fn,'\n',sep='')
    figure.fn = paste(fig.dir,'3dpattern_eig_distr_',clust.name,'.pdf',sep='')
    pdf(figure.fn,useDingbats=F,width=5,height=7)
    print(ggplot(eig.clu.df, aes(x=clust, fill=clust, y=eig)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        guides(fill=FALSE) +
        xlab('') +
        ylab('First eigenvector')
        )
    dev.off()
}


for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]
    abs.gr = abscore.grs[[clust.name]]
    
    cat('[2.4] A/B score (',clust.name,')\n',sep='')
    abs.clust.ov <- findOverlaps(abs.gr, clust.gr, ignore.strand=TRUE)
    abs.gr$clust <- NA
    abs.gr$clust[abs.clust.ov@from] <- clust.gr[abs.clust.ov@to]$global.cluster
    abs.clu.df = data.frame(chr=seqnames(abs.gr), abscore=abs.gr$abscore, clust=abs.gr$clust)
    abs.clu.df = abs.clu.df[!is.na(abs.clu.df$abscore) & !is.na(abs.clu.df$clust),]

    cat('[2.4] Distribution of A/B scores by 3D pattern (',clust.name,')-> ',figure.fn,'\n',sep='')
    figure.fn = paste(fig.dir,'3dpattern_abscore_distr_',clust.name,'.pdf',sep='')
    pdf(figure.fn,useDingbats=F,width=5,height=7)
    print(ggplot(abs.clu.df, aes(x=clust, fill=clust, y=abscore)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        guides(fill=FALSE) +
        xlab('') +
        ylab('A/B score')
        )
    dev.off()
}

##
## 3. Signal distribution in 3D patterns.
##
cat('[3. Signal in 3D patterns]\n')


##
## 3.1. Chromatin marks in 3D patterns.
##
cat('[3.1. Chromatin marks in 3D patterns]\n')

for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]

    ## Preallocate matrix
    mat = list();
    for (chr in levels(seqnames(clust.gr)))
        mat[[chr]] <- matrix(0,ncol=max(clust.gr[seqnames(clust.gr) == chr]$local.cluster),nrow=length(names))
    gwmat <- matrix(0,ncol=length(unique(clust.gr$global.cluster)),nrow=length(names))

    ## Heatmap colors
    hmcols <- colorRampPalette(c("blue","white","red"))(32)

    ## Local intrachromosomal clusters
    normat <- list();
    for (chr in levels(seqnames(clust.gr))) {
        ## Select chromosome
        gc = clust.gr[seqnames(clust.gr) == chr]
        gd = cm.gr[seqnames(cm.gr) == chr]

        ## Overlap Genomic Ranges
        ov = findOverlaps_nw(gc,gd,type='any')
        
        ## Fill matrices
        mat[[chr]] = as.matrix(table(data.frame(dataset=gd[ov@to]$dataset,cluster=gc[ov@from]$local.cluster)))
        clsz = table(gc$local.cluster)
        normat[[chr]] = mat[[chr]]/matrix(clsz,nrow=dim(mat[[chr]])[1],ncol=length(clsz),byrow=TRUE)/rowSums(mat[[chr]])

        ## Intrachromosomal cluster figure
        if (intrachromosomal_chromatin_heatmaps) {
            figure.fn = paste(fig.dir,'heatmap_3dpattern_',clust.name,'_chrommarks_',chr,'.pdf',sep='')
            cat('[3.1] intrachromosomal ',chr,' (',clust.name,') -> ',figure.fn,'\n',sep='')
            pdf(figure.fn,useDingbats=F,width=12,height=8)
            ## Add 'trace = "none"' for presentation figure.
            heatmap.2(normat[[chr]],dendrogram='both',symbreaks=F,labRow=cm.names,scale='row',col=hmcols,margin=c(4,8),
                      symkey=F,key.title="Color key",xlab='3D patterns',
                      main=paste('Jurkat',chr,'3D patterns vs Chromatin marks'))
            dev.off()
        }
    }

    ## Global cluster vs chromatin marks matrix
    clu.cm.ov <- findOverlaps_nw(clust.gr,cm.gr,type='any')
    clu.cm.df <- data.frame(dataset=cm.gr[clu.cm.ov@to]$dataset,cluster=clust.gr[clu.cm.ov@from]$global.cluster)
    clu.cm.gw <- as.matrix(table(clu.cm.df))

    ## Global cluster vs chromatin marks heatmap
    if (interchromosomal_chromatin_heatmap) {
        ## Normalize matrix by rowSums
        clu.cm.mat  <- clu.cm.gw/matrix(clust.bins,nrow=dim(clu.cm.gw)[1],ncol=length(clust.bins),byrow=TRUE)/rowSums(clu.cm.gw)
        ## Create heatmap
        clu.cm.hmap <- heatmap.2(clu.cm.mat,dendrogram='both',symbreaks=F,labRow=cm.names,scale='row',col=hmcols,
                                 margin=c(4,8),symkey=F,key.title="Color key",xlab='3D pattern',
                                 main='Jurkat GW 3D patterns vs Chromatin marks')
        clu.cm.lab  <- colnames(clu.cm.mat)[clu.cm.hmap$colInd]
        ## Save figure
        figure.fn <- paste(fig.dir,'heatmap_3dpattern_',clust.name,'_chrommarks.pdf',sep='')
        cat('[3.1] heatmap: chromatin marks vs 3D patterns (',clust.name,') -> ',figure.fn,'\n',sep='')
        pdf(figure.fn,useDingbats=F,width=12,height=15)
        ## Add 'trace = "none"' for presentation figure
        heatmap.2(clu.cm.mat, dendrogram='both', symbreaks=FALSE, labRow=cm.names,scale='row', col=hmcols,margin=c(4,8),
                  symkey=FALSE, key.title="Color key", xlab='3D pattern', main='Jurkat GW 3D patterns vs Chromatin marks',
                  add.expr = text(x = seq_along(clu.cm.lab), y = 0.1, cex=1.5, srt =0, labels=clu.cm.lab, xpd=TRUE),
                  labCol="")
        dev.off()

        if (interchromosomal_chromatin_heatmaps_light) {
            figure.fn = paste(fig.dir,'heatmap_3dpattern_',clust.name,'_chrommarks_light.pdf',sep='')
            cat('[3.1] heatmap: chromatin marks vs 3D patterns (',clust.name,') -> ',figure.fn,'\n',sep='')
            pdf(figure.fn,useDingbats=F,width=12,height=15)
            ## Add 'trace = "none"' for presentation figure
            heatmap.2(clu.cm.mat, dendrogram='none', symbreaks=FALSE, labRow=cm.names,scale='row',
                      col=hmcols,margin=c(6,10), key=FALSE, xlab='3D pattern', trace="none",
                      main='Jurkat GW 3D patterns vs Chromatin marks',
                      add.expr = text(x = seq_along(clu.cm.lab), y=0.1, cex=1.5, srt=0, labels=clu.cm.lab, xpd=TRUE),
                      labCol="")
            dev.off()
        }
    }

    ## Signal coverage in clusters (barplots)
    clu.cm.cov.df = data.frame(clu=c(), cm=c(), nt=c(), cm.cov=c(), clu.cov=c())
    for (cm.id in unique(cm.gr$dataset)) {
        for (clu in unique(clust.gr$global.cluster)) {
            tmp.df = data.frame(clu    = clu,
                                cm     = cm.names[cm.id],
                                nt     = sum(as.numeric(sum(coverage(intersect(cm.gr[cm.gr$dataset == cm.id], clust.gr[clust.gr$global.cluster == clu])))))
                                )
            tmp.df$cm.cov  <- tmp.df$nt / sum(as.numeric(sum(coverage(cm.gr[cm.gr$dataset == cm.id]))))
            tmp.df$clu.cov <- tmp.df$nt / sum(as.numeric(sum(coverage(clust.gr[clust.gr$global.cluster == clu]))))
            clu.cm.cov.df  <- rbind(clu.cm.cov.df, tmp.df)
        }
    }

    clu.cm.cov.df$clu <- factor(clu.cm.cov.df$clu, c("A1","A2","B1","B2","B3"))

    if (cmarks.type == 'reduced') {
        cm.sname <- unlist(strsplit(as.character(clu.cm.cov.df$cm),'[.]'))
        cm.sname <- cm.sname[seq(2,length(cm.sname),2)]
        clu.cm.cov.df$cm <- factor(cm.sname, c("Pol2","H3K27ac","H3K36me3","H3K4me3","H3K79me3","H3K9me3","H3K27me3","Lamin"))
        clu.cm.cov.df <- clu.cm.cov.df[!is.na(clu.cm.cov.df$cm),]
        fig.w <- 12
        fig.h <- 3
    } else {
        fig.w <- 16
        fig.h <- 12
    }


    figure.fn <- paste(fig.dir,'chromatin_3dpattern_',clust.name,'_clust_coverage.pdf',sep='')
    cat('[3.1] Chromatin signal in 3D patterns (relative cluster coverage, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=fig.w,height=fig.h)
    print(
        ggplot(clu.cm.cov.df, aes(x=cm, y=clu.cov, fill=clu)) +
        geom_bar(stat='identity', position='dodge', colour='black') +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        guides(fill=FALSE) +
        coord_flip() +
        facet_wrap(~clu, ncol=length(global.clust.names)) +
        theme_minimal() +
        xlab('') +
        ylab('Cluster coverage')
    )
    dev.off()

    figure.fn <- paste(fig.dir,'chromatin_3dpattern_',clust.name,'_clust_coverage_flip.pdf',sep='')
    cat('[3.1] Chromatin signal in 3D patterns (relative cluster coverage, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=fig.w,height=fig.h)
    print(
        ggplot(clu.cm.cov.df, aes(x=clu, y=clu.cov, fill=cm)) +
        geom_bar(stat='identity', position='dodge', colour='black') +
        guides(fill=FALSE) +
        facet_wrap(~cm, ncol=length(unique(clu.cm.cov.df$cm))) +
        ##facet_wrap(~cm, ncol=length(unique(clu.cm.cov.df$cm)), scales='free') +
        theme_minimal() +
        xlab('') +
        ylab('Cluster coverage')
    )
    dev.off()


    figure.fn <- paste(fig.dir,'chromatin_3dpattern_',clust.name,'_signal_distribution.pdf',sep='')
    cat('[3.1] Chromatin signal in 3D patterns (distribution of total signal, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=fig.w,height=fig.h)
    print(
        ggplot(clu.cm.cov.df, aes(x=cm, y=cm.cov, fill=clu)) +
        geom_bar(stat='identity', position='dodge', colour='black') +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        coord_flip() +
        facet_wrap(~clu, ncol=length(global.clust.names)) +
        theme_minimal() +
        xlab('') +
        ylab('Signal distribution in clusters [% of total]')
    )
    dev.off()

    figure.fn <- paste(fig.dir,'chromatin_3dpattern_',clust.name,'_nt_coverage.pdf',sep='')
    cat('[3.1] Chromatin signal in 3D patterns (coverage in nt, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=fig.w,height=fig.h)
    print(
        ggplot(clu.cm.cov.df, aes(x=cm, y=nt, fill=clu)) +
        geom_bar(stat='identity', position='dodge', colour='black') +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        coord_flip() +
        facet_wrap(~clu, ncol=length(global.clust.names)) +
        theme_minimal() +
        xlab('') +
        ylab('Signal coverage [bp]')
    )
    dev.off()

    ## A1/A2 enrichment rate
    cccd.a1 = clu.cm.cov.df[clu.cm.cov.df$clu == 'A1',c('cm','clu.cov')]
    cccd.a2 = clu.cm.cov.df[clu.cm.cov.df$clu == 'A2',c('cm','clu.cov')]
    cccd = merge(cccd.a1,cccd.a2,by='cm')
    cccd$enrichment = cccd$clu.cov.x/cccd$clu.cov.y

    figure.fn <- paste(fig.dir,'chromatin_3dpattern_',clust.name,'_clust_coverage_enrichment.pdf',sep='')
    cat('[3.1] Chromatin signal in 3D patterns (A1/A2 cluster coverage enrichment, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=fig.w,height=fig.h)
    print(
        ggplot(cccd, aes(x=cm, y=enrichment, fill=cm)) +
        geom_bar(stat='identity', position='dodge', colour='black') +
        guides(fill=FALSE) +
        coord_flip() +
        theme_minimal() +
        xlab('') +
        ylab('A1/A2 cluster coverage enrichment')
    )
    dev.off()
}

##
## 3.2. HIV in 3D patterns.
##
cat('[3.2. HIV insertion sites in 3D patterns]\n')

hiv.grs = list()
for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]
    hiv.colname = paste('clust.',clust.name,sep='')

    ## Heatmap of HIV insertion sites in 3D patterns
    clu.hiv.ov   <- findOverlaps_nw(clust.gr,hiv.gr,type='any',ignore.strand=T)
    clu.hiv.df   <- data.frame(dataset=hiv.gr[clu.hiv.ov@to]$dataset.name,
                               cluster=clust.gr[clu.hiv.ov@from]$global.cluster)
    clu.hiv.gw   <- as.matrix(table(clu.hiv.df))

    ## Annotate cluster info to HIV integrations.
    tmp.gr                      <- hiv.gr
    tmp.gr$clust                <- NA
    tmp.gr$clust[clu.hiv.ov@to] <- clust.gr[clu.hiv.ov@from]$global.cluster

    hiv.grs[[clust.name]]       <- tmp.gr

    if (interchromosomal_hiv_heatmap) {
        ## Normalize matrix by rowSums
        clu.hiv.mat  <- clu.hiv.gw/matrix(clust.bins,nrow=dim(clu.hiv.gw)[1],ncol=length(clust.bins),byrow=TRUE)/rowSums(clu.hiv.gw)
        ## Create heatmap
        clu.hiv.hmap <- heatmap.2(clu.hiv.mat,dendrogram='column',symbreaks=F,labRow=hiv.names,scale='row',
                                  col=hmcols,margin=c(10,20),symkey=F,key.title="Color key",xlab='3D pattern',
                                  main='Jurkat GW 3D patterns vs HIV insertion sites')
        clu.hiv.lab  <- colnames(clu.hiv.mat)[clu.hiv.hmap$colInd]
        ## Save figure
        figure.fn <- paste(fig.dir,'heatmap_3dpattern_hiv.pdf',sep='')
        cat('[3.2] heatmap: HIV insertion sites vs 3D patterns -> ',figure.fn,'\n',sep='')
        pdf(figure.fn,useDingbats=F,width=16,height=5)
        ## Add 'trace = "none"' for presentation figure
        heatmap.2(clu.hiv.mat, dendrogram='column', symbreaks=FALSE, labRow=hiv.names,scale='row', col=hmcols,
                  margin=c(10,20), key=FALSE, xlab='3D pattern', main='Jurkat GW 3D patterns vs HIV insertion sites',
                  add.expr = text(x = seq_along(clu.hiv.lab), y = 0.1, cex=2, srt =0, labels=clu.hiv.lab,
                                  xpd=TRUE),
                  labCol="")
        dev.off()
    }
}

##
## 4. HIV integration pattern
##
cat('[4. HIV integration pattern]\n');

for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]
    clust.categ.gr = clust.categs.gr[[clust.name]]
    clust.categ = clust.categs[[clust.name]]
    hiv.gr = hiv.grs[[clust.name]]
    
    ## Compute HIV frequencies in cluster/categories.
    clu.cat.hiv.ov    <- findOverlaps_nw(hiv.gr, clust.categ.gr, ignore.strand=TRUE)
    clu.cat.hiv       <- as.data.frame(table(clust.categ.gr[clu.cat.hiv.ov@to]@elementMetadata))
    clu.cat.hiv       <- merge(clust.categ,clu.cat.hiv)
    ## Sort categories and clusters
    clu.cat.hiv$cat   <- factor(clu.cat.hiv$cat, c("actgen","actpro","enh","silgen","int"))
    ##clu.cat.hiv$clust <- factor(clu.cat.hiv$clust, c("Act+","Act-","H3K27me","H3K9me","Lamin"))
    clu.cat.hiv$clust <- factor(clu.cat.hiv$clust, c("A1","A2","B1","B2","B3"))
    clu.cat.hiv       <- clu.cat.hiv[order(clu.cat.hiv$clu, clu.cat.hiv$cat),]

    ##
    ## 4.1. Integration density in Categories.
    ##
    cat('[4.1. Integration density bias in Categories]\n')

    ## Plot spie chart of categories.
    cat.hiv = aggregate(cbind(span, Freq) ~ cat, data=clu.cat.hiv, FUN=sum)

    ## GSpie shapes
    cat.hiv$xmax <- cumsum(cat.hiv$span)
    cat.hiv$xmin <- c(0, cat.hiv[1:(nrow(cat.hiv)-1),]$xmax)
    cat.hiv$ymin <- 0
    cat.hiv$ymax <- sqrt( cat.hiv$Freq/sum(cat.hiv$Freq) / (cat.hiv$span / sum(cat.hiv$span)) )
    ## Plot spie-chart using ggplot2
    figure.fn <- paste(fig.dir,'hiv_bias_category.pdf',sep='')
    cat('[4.1] HIV integration bias in categories -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=7,height=7)
    print(
        ggplot(cat.hiv) +
        geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
        geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
        geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color=categ.colors[cat.hiv$cat]) +
        scale_fill_manual(name='Category',values=categ.colors) +
        coord_polar('x') +
        theme_minimal() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        ggtitle('HIV integration bias in Categories')
    )
    dev.off()


    ##
    ## 4.2. Integration density in 3D patterns.
    ##
    cat('[4.2. Integration density bias in 3D Patterns]\n')

    ## Plot spie chart of categories.
    clu.hiv = aggregate(cbind(span, Freq) ~ clust, data=clu.cat.hiv, FUN=sum)

    ## GSpie shapes
    clu.hiv$xmax <- cumsum(clu.hiv$span)
    clu.hiv$xmin <- c(0, clu.hiv[1:(nrow(clu.hiv)-1),]$xmax)
    clu.hiv$ymin <- 0
    clu.hiv$ymax <- sqrt( clu.hiv$Freq/sum(clu.hiv$Freq) / (clu.hiv$span / sum(clu.hiv$span)) )
    ## Plot spie-chart using ggplot2
    figure.fn <- paste(fig.dir,'hiv_bias_3dpattern_',clust.name,'.pdf',sep='')
    cat('[4.2] HIV integration bias in 3D Patterns (',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=7,height=7)
    print(
        ggplot(clu.hiv) +
        geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
        geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
        geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color='black') +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        coord_polar('x') +
        theme_minimal() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        ggtitle('HIV integration bias in 3D Patterns')
    )
    dev.off()


    ##
    ## 4.3. Local integration density in categories (grouped and normalized by 3D pattern)
    ##
    cat('[4.3. Integration density bias in categories (grouped and normalized by 3D pattern)]\n')

    ## Compute signal/background distributions and spie shapes.
    spie.list = list()
    for (j in unique(clust.gr$global.cluster)) {
        cat.hiv <- clu.cat.hiv[clu.cat.hiv$clust == j,]
        ## GSpie shapes
        cat.hiv$xmax <- cumsum(cat.hiv$span)
        cat.hiv$xmin <- c(0, cat.hiv[1:(nrow(cat.hiv)-1),]$xmax)
        cat.hiv$ymin <- 0
        cat.hiv$ymax <- sqrt( cat.hiv$Freq/sum(cat.hiv$Freq) / (cat.hiv$span / sum(cat.hiv$span)) )
        ## Plot spie-chart using ggplot2
        spie.plt <- ggplot(cat.hiv) +
            geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
            geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
            geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color=categ.colors[cat.hiv$cat]) +
            scale_fill_manual(name='Category',values=categ.colors) +
            coord_polar('x') +
            theme_minimal() +
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                  axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank()) +
            ggtitle(j)
        ## Append to plot list
        spie.list[[j]] = spie.plt
    }

    ## Plot HIV category biases in spie-charts.
    figure.fn <- paste(fig.dir,'hiv_bias_categories_in_3dpatterns_',clust.name,'.pdf',sep='')
    cat('[4.3] HIV local integration bias (category) in 3D patterns (',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=20,height=15)
    print(do.call("grid.arrange", c(spie.list, ncol=floor(sqrt(length(spie.list))))))
    dev.off()


    ##
    ## 4.4. Local integration density in 3D patterns (grouped and normalized by Category)
    ##
    cat('[4.4. Integration density bias in 3D pattern (grouped and normalized by Category)]\n')

    ## Compute signal/background distributions and spie shapes.
    spie.list = list()
    for (i in names(categ)) {
        clu.hiv <- clu.cat.hiv[clu.cat.hiv$cat == i,]
        ## GSpie shapes
        clu.hiv$xmax <- cumsum(clu.hiv$span)
        clu.hiv$xmin <- c(0, clu.hiv[1:(nrow(clu.hiv)-1),]$xmax)
        clu.hiv$ymin <- 0
        clu.hiv$ymax <- sqrt( clu.hiv$Freq/sum(clu.hiv$Freq) / (clu.hiv$span / sum(clu.hiv$span)) )
        ## Plot spie-chart using ggplot2
        spie.plt <- ggplot(clu.hiv) +
            geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
            geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
            geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color='black') +
            scale_fill_manual(name='3D pattern',values=clust.colors) +
            coord_polar('x') +
            theme_minimal() +
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                  axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank()) +
            ggtitle(i)
        ## Append to plot list
        spie.list[[i]] = spie.plt
    }

    ## Plot HIV category biases in spie-charts.
    figure.fn <- paste(fig.dir,'hiv_bias_3dpatterns_',clust.name,'_in_categories.pdf',sep='')
    cat('[4.4] HIV local integration bias (3D pattern, ',clust.name,') in categories -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=20,height=15)
    print(do.call("grid.arrange", c(spie.list, ncol=floor(sqrt(length(spie.list))))))
    dev.off()


    ##
    ## 4.5. HIV integration density (grouped by 3D pattern or category)
    ##
    cat('[4.5. HIV integration density in 3D patterns and categories]\n')

    ## Compute integration density (per Mb) in cluster+category.
    clu.cat.hiv$integ.mb <- clu.cat.hiv$Freq/clu.cat.hiv$span*1e6

    ## Barplot of 3d patterns grouped by category.
    figure.fn <- paste(fig.dir,'hiv_density_category_3dpattern_',clust.name,'.pdf',sep='')
    cat('[4.5] HIV integration density (groups: category, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=8,height=5)
    print(
        ggplot(data=clu.cat.hiv,aes(fill=clust,x=cat,y=integ.mb)) +
        geom_bar(stat='identity',position='dodge',colour='black') +
        scale_fill_manual(name='3D pattern',values=c(clust.colors)) +
        xlab('Category') +
        ylab('HIV integrations per Mbp')
    )
    dev.off()

    ## Barplot of categories grouped by 3d pattern.
    figure.fn <- paste(fig.dir,'hiv_density_3dpattern_',clust.name,'_category.pdf',sep='')
    cat('[4.5] HIV integration density (groups: 3D pattern, ',clust.name,') -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=8,height=5)
    print(
        ggplot(data=clu.cat.hiv,aes(fill=cat,x=clust,y=integ.mb)) +
        geom_bar(stat='identity',position='dodge',colour='black') +
        scale_fill_manual(name='Category',values=c(categ.colors)) +
        xlab('3D patterns') +
        ylab('HIV integrations per Mbp')
    )
    dev.off()

}


##
## 4.6. HIV hotspot detection.
##

hp.count = 5
hp.top.pctl = 5

                                        #for (hp.count in c(10,15,20,50)) {
hp.span = GRanges()
for (c in unique(seqnames(hiv.gr))) {
    hp.chr = sort(hiv.gr[seqnames(hiv.gr) == c], ignore.strand = TRUE)
    if (length(hp.chr) < hp.count) next;
    hp.span = c(hp.span,
                GRanges(seqnames = c,
                        ranges   = IRanges(start = start(hp.chr)[1:(length(hp.chr)-hp.count+1)],end = end(hp.chr)[hp.count:length(hp.chr)])
                        )
                )
}

## Get Top X% dense regions of hp.count.
hp.count.gr = hp.span[order(hp.span@ranges@width)]
hiv.hp.gr   = reduce(hp.count.gr[1:(round(length(hp.span)*hp.top.pctl/100))])
hiv.hp.cnt  = findOverlaps_nw(hiv.hp.gr, hiv.gr, ignore.strand=TRUE)
hiv.hp.cnt  = as.data.frame(table(hiv.hp.cpnt@from))

hiv.hp.gr$hiv.cnt = 0
hiv.hp.gr[as.numeric(hiv.hp.cnt$Var1)]$hiv.cnt = hiv.hp.cnt$Freq
hiv.hp.gr$hiv.dens.kb = hiv.hp.gr$hiv.cnt / width(hiv.hp.gr) * 1000
hiv.hp.gr = hiv.hp.gr[order(hiv.hp.gr$hiv.dens.kb, decreasing=TRUE)]


for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]

    hiv.hp.clu.ov = findOverlaps_nw(hiv.hp.gr, clust.gr, ignore.strand=TRUE)
    hiv.hp.agg = aggregate(clust ~ hp.id, data= data.frame(hp.id = hiv.hp.clu.ov@from, clust=clust.gr[hiv.hp.clu.ov@to]$global.cluster), FUN=MaxTable)
    hiv.hp.gr$clust = NA
    hiv.hp.gr$clust[hiv.hp.agg$hp.id] = hiv.hp.agg$clust

    hp.clu.integ = table(factor(hiv.hp.gr$clust, c('A1','A2','B1','B2','B3')))
    hp.clu.bg    = table(factor(clust.gr$global.cluster, c('A1','A2','B1','B2','B3')))
    hp.clu.df    = data.frame(clu=names(hp.clu.integ),
                              xmin=c(0,cumsum(as.vector(hp.clu.bg))[1:length(hp.clu.bg)-1]),
                              xmax=cumsum(as.vector(hp.clu.bg)),
                              ymin=0,
                              ymax=sqrt((as.vector(hp.clu.integ)/sum(hp.clu.integ))/(as.vector(hp.clu.bg)/sum(hp.clu.bg))),
                              ratio = (as.vector(hp.clu.integ)/sum(hp.clu.integ))/(as.vector(hp.clu.bg)/sum(hp.clu.bg))
                              )

    ## Spie chart of hotspot localization (3D pattern)
    figure.fn <- paste(fig.dir,'hiv_hotspot',hp.count,'_3dpattern_',clust.name,'_spie.pdf',sep='')
    cat('[4.5] 3D pattern (',clust.name,') localization of hotspots (spie chart) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=8,height=5)
    print(
        ggplot(hp.clu.df) +
        geom_hline(yintercept=sqrt(c(0.5,2,5,10)),alpha=0.5,linetype="dotted") +
        geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
        geom_rect(aes(fill=clu,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color='black') +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        coord_polar('x') +
        theme_minimal() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
        ggtitle('3D pattern of HIV hotspots')
    )
    dev.off()

    ## Boxplot of HIV hotspot density (By 3D pattern).
    figure.fn <- paste(fig.dir,'hiv_hotspot',hp.count,'_dens_3dpattern_',clust.name,'_boxplot_A1_A2.pdf',sep='')
    cat('[4.5] hotspot density distribution by 3D pattern (',clust.name,')-> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=5,height=8)
    hiv.hp.gr.a12 = hiv.hp.gr[hiv.hp.gr$clust %in% c('A1','A2')]
    hp.box.df = data.frame(clu = hiv.hp.gr.a12$clust, hiv.dens.kb = hiv.hp.gr.a12$hiv.dens.kb)
    print(
        ggplot(hp.box.df, aes(x = clu, fill = clu, y=hiv.dens.kb)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='3D Pattern',values=clust.colors) +
        theme_minimal() +
        ylab('Hotspot density (integs/kbp)') +
        xlab('3D Pattern')
    )
    dev.off()

    ## Hotspot proximity to Super Enhancer
    hiv.hp.se.ov <- findOverlaps_nw(hiv.hp.gr, sen.gr+5000, ignore.strand=TRUE)
    hiv.hp.gr$sen.close = 'No SuperEnhancer'
    hiv.hp.gr[hiv.hp.se.ov@from]$sen.close = 'Close to SuperEnhancer'

    hiv.hp.sen.df = as.data.frame(table(sen.close=hiv.hp.gr$sen.close))
    hiv.hp.sen.df$xmin = c(0,hiv.hp.sen.df$Freq[1])
    hiv.hp.sen.df$xmax = cumsum(hiv.hp.sen.df$Freq)
    hiv.hp.sen.df$ymin = 0
    hiv.hp.sen.df$ymax = 1

    figure.fn <- paste(fig.dir,'hiv_hotspot',hp.count,'_superenhancer.pdf',sep='')
    cat('[4.5] hotspots and superenhancers -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=8,height=5)
    print(
        ggplot(hiv.hp.sen.df) +
        geom_rect(aes(fill=sen.close, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color='black') +
        coord_polar('x') +
        theme(axis.title.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank()) +
        guides(fill=guide_legend(title='')) +
        ggtitle('HIV hotspots and Super Enhancers (+5kb up/downstream)')
    )
    dev.off()


    ## Distance from SE to closest HP
    sen.clu.ov = findOverlaps_nw(sen.gr, clust.gr, ignore.strand=TRUE)
    sen.clu.idx = data.frame(sen.id = sen.clu.ov@from, clu = clust.gr[sen.clu.ov@to]$global.cluster)
    sen.clu.idx = aggregate(clu ~ sen.id, data = sen.clu.idx, FUN=MaxTable)
    sen.gr$clust = NA
    sen.gr$clust[sen.clu.idx$sen.id] = sen.clu.idx$clu

    sen.gr$dist.to.alg.hp = NA
    sen.hp.dist = distanceToNearest(sen.gr, hiv.hp.gr, ignore.strand=T)
    sen.gr$dist.to.alg.hp[sen.hp.dist@from] = sen.hp.dist@elementMetadata$distance

    ## Distance from SE to closest HIV integ
    sen.gr$dist.to.hiv = NA
    sen.hiv.dist = distanceToNearest(sen.gr, hiv.gr, ignore.strand=T)
    sen.gr$dist.to.hiv[sen.hiv.dist@from] = sen.hiv.dist@elementMetadata$distance

    ##Distance from SE to 3rd nearest HIV.
    hiv.s.gr = sort(hiv.gr, ignore.strand=T)
    hiv3.gr = GRanges()
    for (c in unique(seqnames(hiv.s.gr))) {
        hiv.tmp = hiv.s.gr[seqnames(hiv.s.gr) == c]
        hiv.tmp = hiv.tmp[order(start(hiv.tmp))]
        hiv3.gr = c(hiv3.gr, GRanges(seqnames=c, range=IRanges(start=start(hiv.tmp)[1:(length(hiv.tmp)-3)], end=start(hiv.tmp)[4:length(hiv.tmp)])))
    } 
    sen.hiv3.near = nearest(sen.gr, hiv3.gr, ignore.strand=T)
    sen.hiv3.dist = rep(NA,length(sen.gr))
    sen.gr$hiv3.dist = NA
    for (i in seq(1,length(sen.hiv3.near))) {
        sen.gr$hiv3.dist[i] = width(range(resize(sen.gr[i],width=1,fix="center"),hiv3.gr[sen.hiv3.near[i]]))
    }

    write.table(sen.gr[order(sen.gr$hiv3.dist)],
                file=paste("sen_dist_to_hiv3_",clust.name,".txt",sep=''),
                quote=F,
                sep='\t',
                row.names=F,
                col.names=T)
}

##
## 5. Gene expression in 3D patterns
##
cat('\n[5. Gene expression in 3D patterns]\n')


##
## 5.1 Endogenous expression in 3D patterns
##
cat('[5.1. Endogenous expression in 3D patterns]\n')
exprtss.gr <- resize(genes.gr[genes.gr$expr > 0],1)
noextss.gr <- resize(genes.gr[genes.gr$expr == 0],1)

for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]

    genes.expr.ov  <- findOverlaps_nw(clust.gr, exprtss.gr,ignore.strand=T)
    genes.expr.clu <- data.frame(clust = clust.gr[genes.expr.ov@from]$global.cluster,
                                 expr  = exprtss.gr[genes.expr.ov@to]$expr)

    genes.sil.ov  <- findOverlaps_nw(clust.gr, noextss.gr, ignore.strand=T)
    genes.sil.clu <- data.frame(clust = clust.gr[genes.sil.ov@from]$global.cluster,
                                expr  = noextss.gr[genes.sil.ov@to]$expr)
    ## Sort 3D Pattern factors
    genes.expr.clu$clust <- factor(genes.expr.clu$clust,c("A1","A2","B1","B2","B3"))

    ## Barplot of endogenous expression in 3D patterns.
    figure.fn <- paste(fig.dir,'expr_endog_3dpattern_',clust.name,'.pdf',sep='')
    cat('[5.1] Endogenous expression distribution in 3D patterns (',clust.name,') (boxplot) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=5, height=7)
    print(
        ggplot(genes.expr.clu, aes(x = clust, y = expr, fill = clust)) +
        geom_boxplot(outlier.size=0.3) +
        scale_y_continuous(trans='log10', limits=c(1e-3,1e5), breaks=c(1e-2,1,1e2,1e4), minor_breaks=c(1e-3,1e-1,10,1e3)) +
        scale_fill_manual(name='',values = clust.colors) +
        theme_minimal() +
        xlab('') +
        ylab('Expression of Jurkat active genes [tpm]')
    )
    dev.off()

    ## 3D pattern distribution by expression levels.
    genes.expr.clu$pctl <- ecdf(genes.expr.clu$expr)(genes.expr.clu$expr)
    genes.expr.clu$pctl <- as.factor(ceiling(genes.expr.clu$pctl*10)/10)
    genes.sil.clu$pctl  <- as.factor('Silent')
    genes.all.clu       <- rbind(genes.expr.clu, genes.sil.clu)
    genes.all.clu.cnt   <- as.data.frame(table(genes.all.clu[,c('clust','pctl')]))
    genes.all.clu.cnt$pctl <- factor(genes.all.clu.cnt$pctl, c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,'Silent'))

    ## 3D pattern of genes at TSS by endogenous expression.
    figure.fn <- paste(fig.dir,'gene_3dpatterns_',clust.name,'_by_endogenous_expression.pdf',sep='')
    cat('[5.1] 3D pattern at gene TSS, by gene expression intervals -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    print(
        ggplot(genes.all.clu.cnt, aes(x=pctl, y=Freq, fill=clust)) +
        geom_bar(stat='identity', position='fill') +
        scale_fill_manual(name='', values = clust.colors) +
        xlab(paste('Genes by endogenous expression (Percentile)', sep='')) +
        ylab('3D pattern at TSS') +
        theme_minimal()
    )
    dev.off()
}

##
## 5.2 HIV expression in 3D patterns
##
cat('[5.2. HIV expression in 3D patterns]\n')

for (i in seq(1,length(cluster.list))) {
    clust.name = names(cluster.list)[i]
    clust.gr = clust.grs[[clust.name]]
    clust.categ.gr = clust.categs.gr[[clust.name]]

    ## Compute clusters/categories of BHIVE integrations.
    clu.cat.bhive.ov <- findOverlaps_nw(bhive.gr, clust.categ.gr, ignore.strand=TRUE)
    clu.cat.bhive    <- data.frame(clust = clust.categ.gr[clu.cat.bhive.ov@to]$clust,
                                   cat   = clust.categ.gr[clu.cat.bhive.ov@to]$cat,
                                   expr   = bhive.gr[clu.cat.bhive.ov@from]$expr)
    ## Sort categories and clusters
    clu.cat.bhive$cat   <- factor(clu.cat.bhive$cat, c("actgen","actpro","enh","silgen","int"))
                                        #clu.cat.bhive$clust <- factor(clu.cat.bhive$clust, c("Act+","Act-","H3K27me","H3K9me","Lamin"))
    clu.cat.bhive$clust <- factor(clu.cat.bhive$clust, c("A1","A2","B1","B2","B3"))
    clu.cat.bhive       <- clu.cat.bhive[order(clu.cat.bhive$clu, clu.cat.bhive$cat),]

    ## Expression distribution by clusters.
    p1 <- ggplot(clu.cat.bhive, aes(x = cat, y = expr, fill = cat)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='',values = categ.colors) +
        coord_cartesian(ylim=c(0,3)) +
        theme_minimal() +
        xlab('') +
        ylab('Normalized HIV expression (log10)')

    ## Expression distribution by categories.
    p2 <- ggplot(clu.cat.bhive, aes(x = clust, y = expr, fill = clust)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='',values = clust.colors) +
        coord_cartesian(ylim=c(0,3)) +
        theme_minimal() +
        xlab('') +
        ylab('')

    ## Expression by cluster and category (Grouped by cluster)
    p3 <- ggplot(clu.cat.bhive, aes(x = clust, y = expr, fill = cat)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='Category',values = categ.colors) +
        coord_cartesian(ylim=c(0,3)) +
        theme_minimal() +
        xlab('3D pattern') +
        ylab('Normalized HIV expression (log10)')

    ## Expression by cluster and category (Grouped by category)
    p4 <- ggplot(clu.cat.bhive, aes(x = cat, y = expr, fill = clust)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(name='3D pattern',values = clust.colors) +
        coord_cartesian(ylim=c(0,3)) +
        theme_minimal() +
        xlab('Category') +
        ylab('')

    ## Boxplot figures.
    figure.fn <- paste(fig.dir,'expr_bhive_3dpattern_',clust.name,'.pdf',sep='')
    cat('[5.2] HIV expression in 3D patterns (',clust.name,') (full figure) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=20,height=15)
    print(grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2))
    dev.off()
}

##
## FROM NOW ON, USE DEFAULT CLUSTER
##

def.clust = names(cluster.list)[[1]]
clust.gr = clust.grs[[def.clust]]
hiv.gr = hiv.grs[[def.clust]]
clust.categ = clust.categs[[def.clust]]
clust.categ.gr = clust.categs.gr[[def.clust]]

##
## 6. HIV hotspots
##
cat('[6. HIV hotspots]\n')

## Get HIV integ site list
hp.all        <- data.frame(chrom=hiv.gr@seqnames, locus=hiv.gr@ranges@start)
hp.all$bin    <- floor(hp.all$locus/hp.binsize)
hp.all$bin.hr <- floor(hp.all$locus/hphr.binsize)
## Correct 2x read count for (chrX, chrY)
hp.all$count                                               <- 1
hp.all[hp.all$chrom=='chrX' | hp.all$chrom=='chrY',]$count <- 2
## Append and aggregate counts
hotspots      <- aggregate(count ~ chrom+bin, data=rbind(genome.bins, hp.all[,c('chrom','bin','count')]), FUN=sum)
hotspots.hr   <- aggregate(count ~ chrom+bin.hr, data=rbind(genome.bins.hr, hp.all[,c('chrom','bin.hr','count')]), FUN=sum)
## Compute hotspot percentiles.
hotspots$pctl <- ecdf(hotspots$count)(hotspots$count)
hotspots.hr$pctl <- ecdf(hotspots.hr$count)(hotspots.hr$count)

## Define percentile ranges.
hp.pctl.intv    <- c(0.995,0.99,0.95,0.9,(sum(hotspots$count==0)+1)/sum(hotspots$count>=0),0.0)
hp.pctl.names   <- c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Sparse HIV','No HIV')
hphr.pctl.intv  <- c(0.999,0.995,0.99,0.95,(sum(hotspots.hr$count==0)+1)/sum(hotspots.hr$count>=0),0.0)
hphr.pctl.names   <- c('Top 0.1%','Top 0.1%-0.5%','Top 0.5%-1%','Top 1%-5%','Sparse HIV','No HIV')
## Sort intervals
hp.pctl.ord       <- order(hp.pctl.intv)
hp.pctl.intv.s    <- hp.pctl.intv[hp.pctl.ord]
hp.pctl.names.s   <- hp.pctl.names[hp.pctl.ord]
hphr.pctl.ord     <- order(hphr.pctl.intv)
hphr.pctl.intv.s  <- hphr.pctl.intv[hphr.pctl.ord]
hphr.pctl.names.s <- hphr.pctl.names[hphr.pctl.ord]

## Make GenomicRanges and intersect with cluster/category information.
hp.gr   <- GRanges(seqnames = Rle(hotspots$chrom),
                   ranges   = IRanges(start=hotspots$bin*hp.binsize,end=(hotspots$bin+1)*hp.binsize-1),
                   count    = hotspots$count,
                   pctl     = hotspots$pctl,
                   cluster  = NA,
                   intv     = NA
                   )

hphr.gr <- GRanges(seqnames = Rle(hotspots.hr$chrom),
                   ranges   = IRanges(start=hotspots.hr$bin.hr*hphr.binsize,end=(hotspots.hr$bin.hr+1)*hphr.binsize-1),
                   count    = hotspots.hr$count,
                   pctl     = hotspots.hr$pctl,
                   cluster  = NA,
                   intv     = NA
                   )

##
## 6.0. HIV hotspots, gene expression and 3D pattern list
##

## Compute HIV count per bin
hiv.clu.ov = findOverlaps_nw(clust.gr, hiv.gr, ignore.strand=TRUE)
hiv.clu.df = as.data.frame(table(hiv.clu.ov@from))
hgc.gr = clust.gr
hgc.gr$hiv.cnt = 0
hgc.gr$hiv.cnt[as.numeric(levels(hiv.clu.df$Var1))] = hiv.clu.df$Freq

## Merge with gene info
hiv.gen.clu.ov = findOverlaps(hgc.gr, genes.gr, ignore.strand=TRUE)
hiv.bin.cnt.gen = data.frame(bin=hiv.gen.clu.ov@from, gene=hiv.gen.clu.ov@to, hiv.cnt=hgc.gr$hiv.cnt[hiv.gen.clu.ov@from], clust=hgc.gr$global.cluster[hiv.gen.clu.ov@from])

## Compute max "hiv per bin" and total "hiv cnt" per gene
hgen.max.df = aggregate(hiv.cnt ~ gene, hiv.bin.cnt.gen, FUN=max)
hgen.sum.df = aggregate(hiv.cnt ~ gene, hiv.bin.cnt.gen, FUN=sum)
hgen.clu.df = aggregate(clust ~ gene, hiv.bin.cnt.gen, FUN=head,1)
hgc.df = merge(hgen.max.df, hgen.clu.df, by="gene")
names(hgc.df) = c('gene','hiv.bin.max','clust')
hgc.df = merge(hgc.df, hgen.sum.df, by="gene")
hgc.df = hgc.df[,c('gene','hiv.bin.max','hiv.cnt','clust')]

## Add info to genes.gr
hivngen.gr = genes.gr
hivngen.gr$clust = NA
hivngen.gr$hiv.cnt = 0
hivngen.gr$hiv.bin.max = 0

hivngen.gr$clust[hgc.df$gene] = as.character(hgc.df$clust)
hivngen.gr$hiv.cnt[hgc.df$gene] = hgc.df$hiv.cnt
hivngen.gr$hiv.bin.max[hgc.df$gene] = hgc.df$hiv.bin.max

## Write dataset
hivngen.gr = sort(hivngen.gr)
hivngen.df = data.frame(chr=seqnames(hivngen.gr),
                        beg=start(hivngen.gr),
                        end=end(hivngen.gr),
                        strand=strand(hivngen.gr),
                        gen.name=hivngen.gr$name,
                        active=hivngen.gr$active,
                        expr=hivngen.gr$expr,
                        clust=hivngen.gr$clust,
                        hiv.cnt=hivngen.gr$hiv.cnt,
                        hiv.bin.max=hivngen.gr$hiv.bin.max
                        )
write.table(hivngen.df,
            file      = "genes_expr_cluster_hiv.tsv",
            quote     = FALSE,
            sep       ='\t',
            row.names = FALSE,
            col.names = TRUE
            )

##
## 6.1. Hotspot occupancy
##
cat('[6.1. Hotspot occupancy in 3D patterns]\n')


##
## 6.1.1. Cluster occupancy of HIV hostpots.
##
hp.clu.ov     <- findOverlaps(hp.gr, clust.categ.gr, ignore.strand = TRUE)
hp.clu.ref    <- data.frame(hp.bin = hp.clu.ov@from, clu = clust.categ.gr[hp.clu.ov@to]$clust)
hp.clu.ref    <- aggregate(clu ~ hp.bin, data = hp.clu.ref, FUN=MaxTable)

hp.gr$cluster[hp.clu.ref$hp.bin] <- hp.clu.ref$clu

## Quick cluster assignment to HighResolution hotspots
hphr.clu.ov   <- findOverlaps(hphr.gr, clust.categ.gr, ignore.strand = TRUE)
hphr.gr$cluster[hphr.clu.ov@from] <- clust.categ.gr$clust[hphr.clu.ov@to]


## Classify bins in percentile intervals.
hp.gr$intv    <- factor(hp.pctl.names.s[findInterval(hp.gr$pctl,hp.pctl.intv.s)], hp.pctl.names)
hphr.gr$intv  <- factor(hphr.pctl.names.s[findInterval(hphr.gr$pctl,hphr.pctl.intv.s)], hphr.pctl.names)

## Compute 3D pattern rate for all intervals.
hp.clu.table  <- table(
    data.frame(
        intv  = hp.gr$intv,
        ##        clust = factor(hp.gr$cluster,rev(c("Act+","Act-","H3K27me","H3K9me","Lamin")))
        clust = factor(hp.gr$cluster,rev(c("A1","A2","B1","B2","B3")))
    )
)
hp.clu.rates  <- as.data.frame(hp.clu.table/rowSums(hp.clu.table))

## Structure type: Stacked bar plot
figure.fn <- paste(fig.dir,'hiv_hotspot_3dpatterns.pdf',sep='')
cat('[6.1] Distribution of HIV hotspots in 3D patterns -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(hp.clu.rates, aes(x=intv, y=Freq, fill=clust)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name='', values = clust.colors) +
    xlab(paste('Genomic bins (',round(hp.binsize/1e3),' kb) sorted by HIV density',sep='')) +
    ylab('3D pattern ratio') +
    theme_minimal()
dev.off()

##
## 6.1.2. HIV bins cdf.
##
hp.cdf <- data.frame(x=seq(1,length(hp.gr))*hp.binsize, y=cumsum(sort(hp.gr$count,decreasing=TRUE)))

## HIV aggregation CDF.
figure.fn <- paste(fig.dir,'hiv_hotspot_cdf.pdf',sep='')
cat('[6.1] HIV integration coverage CDF -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(data=hp.cdf,aes(x=x, y=y), color='black') +
    geom_line() +
    theme_minimal() +
    ylab('HIV integrations') +
    xlab('Genome coverage [nt]')
dev.off()


##
## 6.2. HIV expression in hotspots.
##
cat('[6.2. HIV expression in hotspots]\n')

## Cross expression data with hotspot/cluster information.
hp.clu.expr.ov <- findOverlaps_nw(bhive.gr, hp.gr, ignore.strand = TRUE)
hp.clu.expr    <- data.frame(
    intv  = hp.gr$intv[hp.clu.expr.ov@to],
    ##    clust = factor(hp.gr$cluster[hp.clu.expr.ov@to], c("Act+","Act-","H3K27me","H3K9me","Lamin")),
    clust = factor(hp.gr$cluster[hp.clu.expr.ov@to], c("A1","A2","B1","B2","B3")),
    expr  = bhive.gr$expr[hp.clu.expr.ov@from]
)

## Remove No HIV intervals, just in case the integ and expr datasets are not synchronized.
hp.clu.expr.no = hp.clu.expr[hp.clu.expr$intv != "No HIV",]

## Expression of hotspot intervals
figure.fn <- paste(fig.dir,'hiv_expr_hotspot.pdf',sep='')
cat('[6.2] Expression in hotspots -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(hp.clu.expr.no,aes(x=intv,y=expr,fill=intv)) +
    geom_boxplot(outlier.shape=NA, alpha=0.3) +
    xlab('') +
    ylab('HIV normalized expression (log)') +
    ylim(c(0,2)) +
    theme_minimal() +
    guides(fill=guide_legend(title='Hotspots'))
dev.off()

## Expression of hotspot intervals, grouped by 3d pattern (Act+, Act- and Lamin only)
figure.fn <- paste(fig.dir,'hiv_expr_hotspot_3dpattern.pdf',sep='')
cat('[6.2] Expression in hotspots, grouped by 3D Pattern -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=12, height=8)
ggplot(hp.clu.expr[hp.clu.expr$clust %in% c('A1','A2','B3'),],aes(x=clust,y=expr,fill=intv)) +
    geom_boxplot(outlier.shape=NA, alpha=0.3) +
    xlab('') +
    ylab('HIV normalized expression (log)') +
    ylim(c(0,2)) +
    theme_minimal() +
    guides(fill=guide_legend(title='Hotspots'))
dev.off()


## Remove H3K27me (almost no HIV), also remove H3K9me from Top 0.5-1% (poor stats)
hp.clu.tmp <- hp.clu.expr[hp.clu.expr$clust %in% c('A1','A2','B2','B3'),]
hp.clu.tmp <- hp.clu.tmp[!(hp.clu.tmp$clust=='H3K9me' & hp.clu.tmp$intv=='Top 0.5%-1%'),]

## Expression in 3d patterns grouped by interval (The more sparse, the greater the expression difference)
figure.fn <- paste(fig.dir,'hiv_expr_3dpattern_hotspot.pdf',sep='')
cat('[6.2] Expression in 3D patterns, grouped by hotspots -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=12, height=8)
ggplot(hp.clu.tmp, aes(x=intv, y=expr, fill=clust)) +
    geom_boxplot(outlier.shape=NA, alpha=1) +
    scale_fill_manual(values=clust.colors, name='3D Pattern') +
    xlab('') +
    ylab('HIV normalized expression (log)') +
    ylim(c(0,2)) +
    theme_minimal()
dev.off()

##
## 6.3 Hotspots and SuperEnhancers
##
cat('[6.3. Hotspot interaction with SuperEnhancers]\n')

## Hotspot distance to SuperEnhancers (1D)
hp.sen.d  <- distanceToNearest(hp.gr, sen.gr, ignore.strand = TRUE)
hp.sen.df <- data.frame(intv  = hp.gr[hp.sen.d@from]$intv,
                        clust = hp.gr[hp.sen.d@from]$cluster,
                        sen.d = hp.sen.d@elementMetadata$distance/1e6)
## Boxplot 1D distance to SuperEnhancer
figure.fn <- paste(fig.dir,'hiv_hotspot_superenhancer_1d.pdf',sep='')
cat('[6.3] Hotspots 1D distance to closest superenhancer -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=12, height=8)
ggplot(data=hp.sen.df, aes(x=intv, y=sen.d, fill=intv)) +
    geom_boxplot(outlier.shape=NA) +
    scale_y_log10(limits=c(0.01,100)) +
    xlab('Genomic bins (100kb) sorted by HIV density') +
    ylab('Distance to closest Super Enhancer [Mbp]') +
    theme_minimal()
dev.off()


## Write SuperEnhancer ranges
hp.sen.fname = paste(sen.pcs.file,'1',sep='.')
write.table(sen.data[,c('id','chr','beg','end','val','cat')],
            file      = hp.sen.fname,
            quote     = FALSE,
            sep       ='\t',
            row.names = FALSE,
            col.names = TRUE
            )

## Remove scientific notation when printing integer
options(scipen=999)

if (interchromosomal_contact_tests) {
    ## Write intervals
    hp.sen.pcs    <- data.frame(id.1=c(),cnt.1=c(),itv.1=c(),id.2=c(),cnt.2=c(),itv.2=c(),pcs=c())
    for (i in seq(1,length(hp.pctl.names))) {
        cat('[6.3] SuperEnhancers: Computing pcs score for',hp.pctl.names[i],'\n')
        ids.samp     <- which(hp.gr$intv == hp.pctl.names[i])
        ##    ids.samp     <- sample(ids, min(length(ids),hp.samp.size), replace = FALSE)
        hp.sen.samp <-  data.frame(id  = ids.samp,
                                   chr = as.vector(hp.gr[ids.samp]@seqnames),
                                   beg = hp.gr[ids.samp]@ranges@start,
                                   end = hp.gr[ids.samp]@ranges@start + hp.gr[ids.samp]@ranges@width - 1,
                                   val = hp.gr[ids.samp]$count,
                                   itv = hp.gr[ids.samp]$intv
                                   )
        ## Write ranges to compute Hi-C contacts
        hp.fname = paste(hp.file,i,sep='.')
        write.table(hp.sen.samp,
                    file      = hp.fname,
                    quote     = FALSE,
                    sep       ='\t',
                    row.names = FALSE,
                    col.names = TRUE
                    )

        ## Call python script to compute interchromosomal Hi-C contacts between genes.
        hp.sen.outf  = paste(sen.pcs.file,'hp.out',i,sep='.')
        system(paste('python src/pairwise_interaction_score_group.py', hp.sen.fname, hp.fname, cool.file, hic.binsize, hp.sen.outf))

        ## Read result data
        hp.sen.pcs <- rbind(read.table(hp.sen.outf,
                                       sep = '\t',
                                       header = FALSE,
                                       col.names = c('id.1','cnt.1','itv.1','id.2','cnt.2,','itv.2','pcs','reads')
                                       ),
                            hp.sen.pcs)
    }

    ## Compute also PCS for SuperEnhancers with themselves.
    cat('[6.3] SuperEnhancers: Computing self pcs\n')
    hp.sen.outf  = paste(sen.pcs.file,'hp.out.self',sep='.')
    system(paste('python src/pairwise_interaction_score_group.py', hp.sen.fname, hp.sen.fname, cool.file, hic.binsize, hp.sen.outf))

    ## Read result data
    hp.sen.pcs <- rbind(read.table(hp.sen.outf,
                                   sep = '\t',
                                   header = FALSE,
                                   col.names = c('id.1','cnt.1','itv.1','id.2','cnt.2,','itv.2','pcs','reads')
                                   ),
                        hp.sen.pcs)

    ## Compare histograms.
    figure.fn <- paste(fig.dir,'hiv_hotspot_superenhancer.pdf',sep='')
    cat('[6.3] clustering of hotspots with superenhancers -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=12, height=8)
    ggplot(hp.sen.pcs, aes(y=..density..,x=pcs, colour=itv.2)) +
        geom_freqpoly(binwidth=0.00025) +
        scale_y_log10() +
        xlim(c(0,0.0025)) +
        xlab('Interchromosomal pairwise contact score [Hi-C reads/kb^2]') +
        ylab('Density') +
        guides(colour=guide_legend(title='')) +
        ggtitle('Contacts with superenhancers') +
        theme_minimal()
    dev.off()

    ## Barplot for simpler representation.
    hp.sen.pcs.thr  <- 0.001
    hp.sen.pcs.data <- data.frame(name=c(), val=c())
    for (intv in unique(hp.sen.pcs$itv.2)) {
        hp.sen.pcs.data = rbind(hp.sen.pcs.data,
                                data.frame(name= intv,
                                           val = sum(hp.sen.pcs[hp.sen.pcs$itv.2 == intv,'pcs'] >= hp.sen.pcs.thr)/sum(hp.sen.pcs$itv.2 == intv)*100
                                           )
                                )
    }
    ## Sort categories
    hp.sen.pcs.data$name <- factor(hp.sen.pcs.data$name, c('SuperEnhancer','Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Sparse HIV','No HIV'))
    ## Thresholded IPCS barplot
    figure.fn <- paste(fig.dir,'hiv_hotspot_superenhancer_barplot.pdf',sep='')
    cat('[6.3] clustering of hotspots with superenhancers (barplot, thr =', hp.sen.pcs.thr,' reads/kb^2) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(hp.sen.pcs.data, aes(y=val,x=name, fill=name)) +
        geom_bar(stat='identity',colour='black',alpha=0.6) +
        scale_fill_manual(values=c('black',hue_pal()(6)), name='') +
        ylab(paste('Interactions with IPCS >=',hp.sen.pcs.thr,'reads/kb^2 [%]')) +
        xlab('') +
        ggtitle('Contacts with SuperEnhancers') +
        theme_minimal()
    dev.off()

}

##
## 6.4. Chromatin composition in Hotspots
##
cat('[6.4. Chromatin composition in hotspots]\n')

## Compute intersect between hostpot and chromatin marks.
hp.cm.df = data.frame(intv=c(), cm=c(), nt=c(), intv.cov=c(), cm.cov=c(), enrich = c(), cm.bins = c(), bins.cov = c())
for (cm.id in unique(cm.gr$dataset)) {
    for (hpctl in hphr.pctl.names) {
        tmp.df = data.frame(intv = hpctl,
                            cm   = cm.names[cm.id],
                            nt   = sum(as.numeric(sum(coverage(intersect(cm.gr[cm.gr$dataset == cm.id], hphr.gr[hphr.gr$intv == hpctl])))))
                            )
        cm.cov          <- sum(as.numeric(sum(coverage(cm.gr[cm.gr$dataset == cm.id]))))
        intv.cov        <- sum(as.numeric(sum(coverage(hphr.gr[hphr.gr$intv == hpctl]))))
        tmp.df$cm.cov   <- tmp.df$nt / cm.cov
        tmp.df$intv.cov <- tmp.df$nt / intv.cov
        tmp.df$enrich   <- (tmp.df$nt / cm.cov) / (intv.cov / sum(as.numeric(sum(coverage(hphr.gr)))))
        tmp.df$cm.bins  <- sum(countOverlaps(hphr.gr[hphr.gr$intv == hpctl], cm.gr[cm.gr$dataset == cm.id]) > 0)
        tmp.df$bins.cov <- tmp.df$cm.bins/length(hphr.gr[hphr.gr$intv == hpctl]) * 100
        hp.cm.df        <- rbind(hp.cm.df, tmp.df)
    }
}

## Barplot hotspot relative coverage
figure.fn <- paste(fig.dir,'hiv_hotspot_chromatin_intv_coverage.pdf',sep='')
cat('[6.4] Chromatin composition in hotspots (relative hotspot coverage) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=16,height=12)
ggplot(hp.cm.df, aes(x=cm, y=intv.cov, fill=intv)) +
    geom_bar(stat='identity', position='dodge', colour='black') +
    coord_flip() +
    guides(fill=guide_legend(title='')) +
    facet_wrap(~intv, ncol=length(unique(hp.cm.df$intv))) +
    theme_minimal() +
    xlab('') +
    ylab('Hotspot coverage')
dev.off()

## Barplot enrichment over expected (if homogeneously distributed)
figure.fn <- paste(fig.dir,'hiv_hotspot_chromatin_enrichment.pdf',sep='')
cat('[6.4] Chromatin composition in hotspots (enrichment over expected) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=16,height=12)
ggplot(hp.cm.df, aes(x=cm, y=enrich, fill=intv)) +
    geom_bar(stat='identity', position='dodge', colour='black') +
    coord_flip() +
    scale_y_log10() +
    guides(fill=guide_legend(title='')) +
    facet_wrap(~intv, ncol=length(unique(hp.cm.df$intv))) +
    theme_minimal() +
    xlab('') +
    ylab('Signal enrichment over homogeneous distribution')
dev.off()

## Bins with signal (rather than nt coverage)
figure.fn <- paste(fig.dir,'hiv_hotspot_chromatin_bin_coverage.pdf',sep='')
cat('[6.4] Chromatin composition in hotspots (hotspots bin coverage) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=16,height=12)
ggplot(hp.cm.df, aes(x=cm, y=bins.cov, fill=intv)) +
    geom_bar(stat='identity', position='dodge', colour='black') +
    coord_flip() +
    guides(fill=guide_legend(title='')) +
    facet_wrap(~intv, ncol=length(unique(hp.cm.df$intv))) +
    theme_minimal() +
    xlab('') +
    ylab('Hotspot bins with enriched signal [%]')
dev.off()

##
## 6.5. Hotspots controls (PolII, H3K36me3)
##
cat('[6.5. Hotspot controls]\n')

if (interchromosomal_contact_tests) {
    ##
    ## 6.5.1. H3K36me3
    ##
    cat('[6.5.1. Hotspot control: H3K36me3]\n')

    ## Read signal
    hp.ctrl.data <- read.table(hp.ctrl.hk36, col.names=c('chr','beg','end','name','value','strand','V1','V2','signif'))
    hp.ctrl.gr   <- GRanges(seqnames = Rle(hp.ctrl.data$chr),
                            ranges   = IRanges(start=hp.ctrl.data$beg, end=hp.ctrl.data$end),
                            value    = hp.ctrl.data$value,
                            signif   = hp.ctrl.data$signif)
    ## Prepare genome bins.
    hp.bin.gr    <- unlist(tile(genome.gr, width=hp.binsize))
    ## Find overlap between ChIP signal and genome bins.
    hp.ctrl.bin.ov <- findOverlaps_nw(hp.ctrl.gr, hp.bin.gr, ignore.strand = TRUE)
    hp.ctrl.cnt    <- aggregate(value ~ bin.id,
                                data = data.frame(bin.id = as.integer(hp.ctrl.bin.ov@to),
                                                  value = hp.ctrl.gr$value[hp.ctrl.bin.ov@from]),
                                FUN = sum)
    ## Assign signal count
    hp.bin.gr$cnt                     <- 0
    hp.bin.gr$cnt[hp.ctrl.cnt$bin.id] <- hp.ctrl.cnt$value
    ## Sort by signal count
    hp.bin.gr      <- hp.bin.gr[order(hp.bin.gr$cnt, decreasing=T)]
    ## Compute hotspot percentiles.
    hp.bin.gr$pctl <- ecdf(hp.bin.gr$cnt)(hp.bin.gr$cnt)

    ## Define percentile ranges.
    ctrl.pctl.intv    <- c(0.995,0.99,0.95,0.9,(sum(hp.bin.gr$cnt==0)+1)/sum(hp.bin.gr$cnt>=0),0.0)
    ctrl.pctl.names   <- c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Weak peaks','No signal')
    ctrl.pctl.colors  <- hue_pal()(length(ctrl.pctl.names))
    ## Sort intervals
    ctrl.pctl.ord     <- order(ctrl.pctl.intv)
    ctrl.pctl.intv.s  <- ctrl.pctl.intv[ctrl.pctl.ord]
    ctrl.pctl.names.s <- ctrl.pctl.names[ctrl.pctl.ord]

    ## Find intervals
    hp.bin.gr$intv <- factor(ctrl.pctl.names.s[findInterval(hp.bin.gr$pctl,ctrl.pctl.intv.s)], ctrl.pctl.names)
    hp.hk36.pcs    <- data.frame(id.1=c(),cnt.1=c(),itv.1=c(),id.2=c(),cnt.2=c(),itv.2=c(),pcs=c())

    ## Remove scientific notation when printing integer
    options(scipen=999)
    ## Compute contact scores for each interval. (use random sampling)
    for (i in seq(1,length(ctrl.pctl.names))) {
        cat('[6.5.1] H3K36me3 control: Computing pcs score for',ctrl.pctl.names[i],'\n')
        ids          <- which(hp.bin.gr$intv == ctrl.pctl.names[i])
        ids.samp     <- sample(ids, min(length(ids),hp.samp.size), replace = FALSE)
        hp.ctrl.samp <-  data.frame(id  = ids.samp,
                                    chr = as.vector(hp.bin.gr[ids.samp]@seqnames),
                                    beg = hp.bin.gr[ids.samp]@ranges@start,
                                    end = hp.bin.gr[ids.samp]@ranges@start + hp.bin.gr[ids.samp]@ranges@width - 1,
                                    val = hp.bin.gr[ids.samp]$cnt,
                                    itv = hp.bin.gr[ids.samp]$intv
                                    )
        ## Write ranges to compute Hi-C contacts
        hp.ctrl.fname = paste(hp.samp.file,i,sep='.')
        write.table(hp.ctrl.samp,
                    file      = hp.ctrl.fname,
                    quote     = FALSE,
                    sep       ='\t',
                    row.names = FALSE,
                    col.names = TRUE
                    )

        ## Call python script to compute interchromosomal Hi-C contacts between genes.
        system(paste('python src/pairwise_interaction_score.py', hp.ctrl.fname, cool.file, hic.binsize))

        ## Read result data
        hp.hk36.pcs <- rbind(read.table(paste(hp.ctrl.fname,'pcs',sep='.'),
                                        sep = '\t',
                                        header = FALSE,
                                        col.names = c('id.1','cnt.1','itv.1','id.2','cnt.2,','itv.2','pcs','reads')
                                        ),
                             hp.hk36.pcs)
    }

    ## Compare histograms.
    figure.fn <- paste(fig.dir,'hiv_hotspot_control_H3K36me3.pdf',sep='')
    cat('[6.5.1] Hotspot control: clustering of H3K36me3 dense regions -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=12, height=8)
    ggplot(hp.hk36.pcs, aes(y=..density..,x=pcs, colour=itv.1)) +
        geom_freqpoly(binwidth=0.00025) +
        scale_y_log10() +
        xlim(c(0,0.002)) +
        xlab('Interchromosomal pairwise contact score [Hi-C reads/kb^2]') +
        ylab('Density') +
        guides(colour=guide_legend(title='')) +
        ggtitle('Control: H3K36me3 clustering') +
        theme_minimal()
    dev.off()

    ## Barplot for simpler representation.
    hp.hk36.pcs.thr  <- 0.0005
    hp.hk36.pcs.data <- data.frame(name=c(), val=c())
    for (intv in unique(hp.hk36.pcs$itv.2)) {
        hp.hk36.pcs.data = rbind(hp.hk36.pcs.data,
                                 data.frame(name= intv,
                                            val = sum(hp.hk36.pcs[hp.hk36.pcs$itv.2 == intv,'pcs'] >= hp.hk36.pcs.thr)/sum(hp.hk36.pcs$itv.2 == intv)*100
                                            )
                                 )
    }
    ## Sort categories
    hp.hk36.pcs.data$name <- factor(hp.hk36.pcs.data$name, c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Weak peaks','No signal'))
    ## Thresholded IPCS barplot
    figure.fn <- paste(fig.dir,'hiv_hotspot_control_H3K36me3_barplot.pdf',sep='')
    cat('[6.5.1] Hotspot control: clustering of H3K36me3 dense regions (barplot, thr =', hp.hk36.pcs.thr,' reads/kb^2) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(hp.hk36.pcs.data, aes(y=val,x=name, fill=name)) +
        geom_bar(stat='identity',colour='black',alpha=0.6) +
        ylab(paste('Interactions with IPCS >=',hp.hk36.pcs.thr,'reads/kb^2 [%]')) +
        xlab('') +
        guides(fill=guide_legend(title='')) +
        ggtitle('Contact score of H3K36me3 sites') +
        theme_minimal()
    dev.off()


    ##
    ## 6.5.2. PolII
    ##
    cat('[6.5.2. Hotspot control: PolII]\n')

    ## Read signal
    hp.ctrl.data <- read.table(hp.ctrl.pol2, col.names=c('chr','beg','end','name','score','strand','value','pvalue','qvalue','peak'))
    hp.ctrl.gr   <- GRanges(seqnames = Rle(hp.ctrl.data$chr),
                            ranges   = IRanges(start=hp.ctrl.data$beg, end=hp.ctrl.data$end),
                            value    = hp.ctrl.data$value,
                            signif   = hp.ctrl.data$qvalue)
    ## Find overlap between ChIP signal and genome bins.
    hp.ctrl.bin.ov <- findOverlaps_nw(hp.ctrl.gr, hp.bin.gr, ignore.strand = TRUE)
    hp.ctrl.cnt    <- aggregate(value ~ bin.id,
                                data = data.frame(bin.id = as.integer(hp.ctrl.bin.ov@to),
                                                  value = hp.ctrl.gr$value[hp.ctrl.bin.ov@from]),
                                FUN = sum)
    ## Assign signal count
    hp.bin.gr$cnt                     <- 0
    hp.bin.gr$cnt[hp.ctrl.cnt$bin.id] <- hp.ctrl.cnt$value
    ## Sort by signal count
    hp.bin.gr      <- hp.bin.gr[order(hp.bin.gr$cnt, decreasing=T)]
    ## Compute hotspot percentiles.
    hp.bin.gr$pctl <- ecdf(hp.bin.gr$cnt)(hp.bin.gr$cnt)

    ## Define percentile ranges.
    ctrl.pctl.intv    <- c(0.995,0.99,0.95,0.9,(sum(hp.bin.gr$cnt==0)+1)/sum(hp.bin.gr$cnt>=0),0.0)
    ctrl.pctl.names   <- c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Weak peaks','No signal')
    ctrl.pctl.colors  <- hue_pal()(length(ctrl.pctl.colors))
    ## Sort intervals
    ctrl.pctl.ord     <- order(ctrl.pctl.intv)
    ctrl.pctl.intv.s  <- ctrl.pctl.intv[ctrl.pctl.ord]
    ctrl.pctl.names.s <- ctrl.pctl.names[ctrl.pctl.ord]

    ## Find intervals
    hp.bin.gr$intv <- factor(ctrl.pctl.names.s[findInterval(hp.bin.gr$pctl,ctrl.pctl.intv.s)], ctrl.pctl.names)
    hp.pol2.pcs    <- data.frame(id.1=c(),cnt.1=c(),itv.1=c(),id.2=c(),cnt.2=c(),itv.2=c(),pcs=c())

    ## Remove scientific notation when printing integer
    options(scipen=999)
    ## Compute contact scores for each interval. (use random sampling)
    for (i in seq(1,length(ctrl.pctl.names))) {
        cat('[6.5.2] PolII control: Computing pcs score for',ctrl.pctl.names[i],'\n')
        ids          <- which(hp.bin.gr$intv == ctrl.pctl.names[i])
        ids.samp     <- sample(ids, min(length(ids),hp.samp.size), replace = FALSE)
        hp.ctrl.samp <-  data.frame(id  = ids.samp,
                                    chr = as.vector(hp.bin.gr[ids.samp]@seqnames),
                                    beg = hp.bin.gr[ids.samp]@ranges@start,
                                    end = hp.bin.gr[ids.samp]@ranges@start + hp.bin.gr[ids.samp]@ranges@width - 1,
                                    val = hp.bin.gr[ids.samp]$cnt,
                                    itv = hp.bin.gr[ids.samp]$intv
                                    )
        ## Write ranges to compute Hi-C contacts
        hp.ctrl.fname = paste(hp.samp.file,i,sep='.')
        write.table(hp.ctrl.samp,
                    file      = hp.ctrl.fname,
                    quote     = FALSE,
                    sep       ='\t',
                    row.names = FALSE,
                    col.names = TRUE
                    )

        ## Call python script to compute interchromosomal Hi-C contacts between genes.
        system(paste('python src/pairwise_interaction_score.py', hp.ctrl.fname, cool.file, hic.binsize))

        ## Read result data
        hp.pol2.pcs <- rbind(read.table(paste(hp.ctrl.fname,'pcs',sep='.'),
                                        sep = '\t',
                                        header = FALSE,
                                        col.names = c('id.1','cnt.1','itv.1','id.2','cnt.2,','itv.2','pcs','reads')
                                        ),
                             hp.pol2.pcs)
    }

    ## Compare histograms.
    figure.fn <- paste(fig.dir,'hiv_hotspot_control_PolII.pdf',sep='')
    cat('[6.5.2] Hotspot control: clustering of PolII dense regions -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=12, height=8)
    ggplot(hp.pol2.pcs, aes(y=..density..,x=pcs, colour=itv.1)) +
        geom_freqpoly(binwidth=0.00025) +
        scale_y_log10() +
        xlim(c(0,0.002)) +
        xlab('Interchromosomal pairwise contact score [Hi-C reads/kb^2]') +
        ylab('Density') +
        guides(colour=guide_legend(title='')) +
        ggtitle('Control: PolII clustering') +
        theme_minimal()
    dev.off()

    ## Barplot for simpler representation.
    hp.pol2.pcs.thr  <- 0.0005
    hp.pol2.pcs.data <- data.frame(name=c(), val=c())
    for (intv in unique(hp.pol2.pcs$itv.2)) {
        hp.pol2.pcs.data = rbind(hp.pol2.pcs.data,
                                 data.frame(name= intv,
                                            val = sum(hp.pol2.pcs[hp.pol2.pcs$itv.2 == intv,'pcs'] >= hp.pol2.pcs.thr)/sum(hp.pol2.pcs$itv.2 == intv)*100
                                            )
                                 )
    }
    ## Sort categories
    hp.pol2.pcs.data$name <- factor(hp.pol2.pcs.data$name, c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Weak peaks','No signal'))
    ## Thresholded IPCS barplot
    figure.fn <- paste(fig.dir,'hiv_hotspot_control_PolII_barplot.pdf',sep='')
    cat('[6.5.2] Hotspot control: clustering of PolII dense regions (barplot, thr =', hp.pol2.pcs.thr,' reads/kb^2) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(hp.pol2.pcs.data, aes(y=val,x=name, fill=name)) +
        geom_bar(stat='identity',colour='black',alpha=0.6) +
        ylab(paste('Interactions with IPCS >=',hp.pol2.pcs.thr,'reads/kb^2 [%]')) +
        xlab('') +
        guides(fill=guide_legend(title='')) +
        ggtitle('Contact score of PolII sites') +
        theme_minimal()
    dev.off()


}

##
## 7. Gene ontology in 3D patterns
##
                                        #genes.clu.ov <- findOverlaps(genes.gr, clust.gr, ignore.strand=TRUE)
                                        #genes.clu    <- data.frame(gene = genes.gr[genes.clu.ov@from]$name, clust = clust.gr[genes.clu.ov@to]$global.cluster)
                                        #genes.clu    <- aggregate(clust ~ gene, genes.clu, FUN=MaxTable)



##
## 8. RIG Genes in Jurkat
##
cat('[8. HIV-dense genes (HD genes)]\n')

##
## 8.1. Clustering of HIV-dense active genes
##
cat('[8.1. Clustering of HD genes]\n')
hiv.gen.ov   <- findOverlaps_nw(hiv.gr, genes.gr, ignore.strand=T)
hiv.gen.cnt  <- aggregate(cnt ~ gene.id, data=data.frame(gene.id=hiv.gen.ov@to, cnt=1), FUN=length)
hiv.gen      <- data.frame(name = genes.gr$name,
                           nt   = genes.gr@ranges@width,
                           cnt  = 0,
                           chr  = as.vector(seqnames(genes.gr)),
                           beg  = genes.gr@ranges@start,
                           end  = genes.gr@ranges@start + genes.gr@ranges@width - 1,
                           act  = genes.gr$active
                           )

## Assign integration count
hiv.gen[hiv.gen.cnt$gene.id,'cnt'] <- hiv.gen.cnt$cnt
hiv.gen$dens.kb                    <- hiv.gen$cnt/hiv.gen$nt*1e3
hiv.gen                            <- subset(hiv.gen, grepl("^chr[0-9XY][0-9]*$",chr))
hiv.gen                            <- hiv.gen[order(hiv.gen$dens.kb, decreasing=T),]

## Write list of HIV-dense genes (+backgroud) to compute pairwise interaction score from Hi-C data.
hiv.ag.rig     <- head(hiv.gen[hiv.gen$nt >= hiv.gen.len & hiv.gen$act, c('name','chr','beg','end','dens.kb','act')], hiv.gen.samp)
hiv.sg.rig     <- head(hiv.gen[hiv.gen$nt >= hiv.gen.len & !hiv.gen$act, c('name','chr','beg','end','dens.kb','act')], hiv.gen.samp)

## Now take a sample with the same number of genes per chromosome.
hiv.ag.rig.chr <- table(hiv.ag.rig$chr)
hiv.sg.rig.chr <- table(hiv.ag.rig$chr)
hiv.ag.bg      <- data.frame(name=c(),chr=c(),beg=c(),end=c(),dens.kb=c())
hiv.sg.bg      <- data.frame(name=c(),chr=c(),beg=c(),end=c(),dens.kb=c())
for (i in seq(1,length(hiv.ag.rig.chr))) {
    hiv.ag.bg  <- rbind(hiv.ag.bg,
                        hiv.gen[sample(which(hiv.gen$cnt == 0 & hiv.gen$act & hiv.gen$chr == names(hiv.ag.rig.chr)[i]),
                                       hiv.ag.rig.chr[i], replace=FALSE),
                                c('name','chr','beg','end','dens.kb','act')]
                        )
    hiv.sg.bg  <- rbind(hiv.sg.bg,
                        hiv.gen[sample(which(hiv.gen$cnt == 0 & !hiv.gen$act & hiv.gen$chr == names(hiv.sg.rig.chr)[i]),
                                       hiv.sg.rig.chr[i], replace=FALSE),
                                c('name','chr','beg','end','dens.kb','act')]
                        )

}

if (interchromosomal_contact_tests) {
    ## Write gene file for processing.
    write.table(rbind(hiv.ag.rig, hiv.ag.bg, hiv.sg.rig, hiv.sg.bg),
                file      = hiv.gen.file,
                quote     = FALSE,
                sep       ='\t',
                row.names = FALSE,
                col.names = TRUE
                )

    ## Call python script to compute interchromosomal Hi-C contacts between genes.
    system(paste('python src/pairwise_interaction_score.py',hiv.gen.file, cool.file, hic.binsize))

    hiv.gen.int <- read.table(paste(hiv.gen.file,'pcs',sep='.'),
                              col.names = c('name.1','dens.1','act.1','name.2','dens.2,','act.2','pcs','reads')
                              )

    cat.names       <- list('No HIV-No HIV','No HIV-RIG','RIG-RIG')
    hiv.gen.int$cat <- unlist(cat.names[as.numeric(hiv.gen.int$dens.1 > 0) + as.numeric(hiv.gen.int$dens.2 > 0)+1])
    hiv.gen.int$cat <- factor(hiv.gen.int$cat, c('RIG-RIG','No HIV-RIG','No HIV-No HIV'))
    act.names       <- list('Silent-Silent','Active-Silent','Active-Active')
    hiv.gen.int$act <- unlist(act.names[as.numeric(hiv.gen.int$act.1) + as.numeric(hiv.gen.int$act.2)+1])
    hiv.gen.int$act <- factor(hiv.gen.int$act, c('Active-Active','Active-Silent','Silent-Silent'))


    ## PCS Density plot
    figure.fn <- paste(fig.dir,'hiv_density_gene_expression_cluster.pdf',sep='')
    cat('[8.1] Clustering of HIV-dense genes by activity and density -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=12, height=8)
    ggplot(hiv.gen.int, aes(y=..density..,x=pcs, colour=act, linetype=cat)) +
        geom_freqpoly(binwidth=0.001) +
        scale_y_log10() +
        xlim(c(0,0.01)) +
        xlab('Interchromosomal pairwise contact score [Hi-C reads / kb^2]') +
        ylab('density') +
        guides(color=guide_legend(title='Gene activity'), linetype=guide_legend('Gene HIV-density')) +
        theme_minimal()
    dev.off()

    ## Compute values for gene activity/density PCS barplot
    hiv.gen.int.pcs.thr  <- 0.0025
    hiv.gen.int.pcs.data <- data.frame(act=c(), cat=c(), val=c())
    for (cat in unique(hiv.gen.int$cat)) {
        for (act in unique(hiv.gen.int$act)) {
            hiv.gen.int.pcs.data = rbind(hiv.gen.int.pcs.data,
                                         data.frame(act = act,
                                                    cat = cat,
                                                    val = sum(hiv.gen.int[hiv.gen.int$cat == cat &
                                                                          hiv.gen.int$act == act,'pcs'] >= hiv.gen.int.pcs.thr)/sum(hiv.gen.int$cat == cat & hiv.gen.int$act == act)*100
                                                    )
                                         )
        }
    }
    hiv.gen.int.pcs.data$act <- factor(hiv.gen.int.pcs.data$act, c('Active-Active','Active-Silent','Silent-Silent'))

    ## Barplot
    figure.fn <- paste(fig.dir,'hiv_density_gene_expression_cluster_barplot.pdf',sep='')
    cat('[8.1] Clustering of HIV-dense genes by activity and density (barplot, thr = ',hiv.gen.int.pcs.thr,' reads/kb^2) -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(hiv.gen.int.pcs.data, aes(y=val,x=cat, fill=act)) +
        geom_bar(stat='identity',position='dodge',colour='black', alpha=0.6) +
        ylab(paste('Interactions with IPCS >=',hiv.gen.int.pcs.thr,'reads/kb^2 [%]')) +
        xlab('') +
        guides(fill=guide_legend(title='')) +
        ggtitle('Contacts between endogenous genes') +
        theme_minimal()
    dev.off()
}

##
## 8.1.2. Clustering of HD gene hotspots (5 kb bins inside genes).
##

## Merge cluster data (5kb on mappable regions) with HIV integration sites.
hiv.bin.ov  <- findOverlaps_nw(clust.gr, hiv.gr, ignore.strand=TRUE)
hiv.bin.cnt <- as.data.frame(table(hiv.bin.ov@from))
hivbin.gr   <- clust.gr

hivbin.gr$hiv.cnt <- 0
hivbin.gr$hiv.cnt[as.numeric(as.character(hiv.bin.cnt[,1]))] = as.numeric(as.character(hiv.bin.cnt[,2]))

genes.sen.ov <- findOverlaps_nw(resize(genes.gr, width(genes.gr)+5000, fix='end'), sen.gr)
genes.gr$sen <- FALSE
genes.gr$sen[unique(genes.sen.ov@from)] <- TRUE

## Overlap between gene bodies and all cluster bins (HIV and No-HIV).
hiv.bin.gen.ov <- findOverlaps_nw(hivbin.gr, genes.gr, ignore.strand=TRUE)

hivbin.gr$gene <- FALSE
hivbin.gr$gene[unique(hiv.bin.gen.ov@from)] <- TRUE

hivbin.gr$active <- NA
hivbin.gr$active[hiv.bin.gen.ov@from] <- genes.gr$active[hiv.bin.gen.ov@to]

hivbin.gr$sen <- FALSE
hivbin.gr$sen[hiv.bin.gen.ov@from] <- genes.gr$sen[hiv.bin.gen.ov@to]

hivbingen.gr <- hivbin.gr[unique(hiv.bin.gen.ov@from)]
hivbingen.gr$gene.expr = 'Silent'
hivbingen.gr$gene.expr[hivbingen.gr$active] = 'Active'

hivbingen.gr$sen.close = 'No SE'
hivbingen.gr$sen.close[hivbingen.gr$sen] = 'SE'

if (interchromosomal_contact_tests) {
    hivbingen.df <- data.frame(chr = seqnames(hivbingen.gr),
                               beg = start(hivbingen.gr),
                               end = end(hivbingen.gr),
                               hiv = hivbingen.gr$hiv.cnt,
                               act = hivbingen.gr$gene.expr)

    hiv.gen.bin.file <- paste(hiv.gen.file,'bins',sep='.')
    hiv.gen.bin.outf <- paste(hiv.gen.file,'bins.out',sep='.')
    write.table(hivbingen.df,
                file      = hiv.gen.bin.file,
                quote     = FALSE,
                sep       ='\t',
                row.names = FALSE,
                col.names = TRUE
                )

    system(paste('python src/pairwise_interaction_bin_groups.py',hiv.gen.bin.file, cool.file, hic.binsize, hiv.gen.bin.outf))

    hiv.gen.bin.pcs <- read.table(hiv.gen.bin.outf, sep='\t', col.names=c('chr.a','act.a','hiv.a','span.a','chr.b','act.b','hiv.b','span.b','pcs','hic.cnt'))

    hgbpcs <- hiv.gen.bin.pcs

    hgbpcs$act <- NA
    hgbpcs$act[hgbpcs$act.a == hgbpcs$act.b & hgbpcs$act.a == 'Active'] = 'Active-Active'
    hgbpcs$act[hgbpcs$act.a == hgbpcs$act.b & hgbpcs$act.a == 'Silent'] = 'Silent-Silent'
    hgbpcs$act[hgbpcs$act.a != hgbpcs$act.b] = 'Active-Silent'

    hgbpcs$hiv <- NA
    hgbpcs$hiv[hgbpcs$hiv.a == hgbpcs$hiv.b & hgbpcs$hiv.a == 'HIV'] = 'HIV-HIV'
    hgbpcs$hiv[hgbpcs$hiv.a == hgbpcs$hiv.b & hgbpcs$hiv.a == 'No HIV'] = 'No HIV-No HIV'
    hgbpcs$hiv[hgbpcs$hiv.a != hgbpcs$hiv.b] = 'HIV-No HIV'

    hgbpcs.rev <- hgbpcs[,c('chr.b','act.a','hiv.a','span.a','chr.a','act.b','hiv.b','span.b','pcs','hic.cnt','act','hiv')]
    names(hgbpcs.rev) <- c('chr.a','act.a','hiv.a','span.a','chr.b','act.b','hiv.b','span.b','pcs','hic.cnt','act','hiv')
    hgbpcs.dup <- rbind(hgbpcs, hgbpcs.rev)

    ## Boxplot
    figure.fn <- paste(fig.dir,'pcs_gene_bins_by_expr_and_hiv.pdf',sep='')
    cat('[8.1.2] Clustering of bins in act/silent genes with/out HIV -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(data=hgbpcs, aes(x=act, y=pcs, fill=hiv)) +
        geom_boxplot(outlier.shape=NA) +
        ylim(c(0,4e-4)) +
        ##    scale_y_log10() +
        theme_minimal() +
        guides(fill=guide_legend(title='')) +
        xlab('Gene activity') +
        ylab('Interchromosomal Hi-C contact density [Hi-C reads/kb^2]')
    dev.off()

    ## Barplot
    hgbpcs.thr <- 2e-4
    hgbpcs.tot = as.data.frame(table(hgbpcs[,c('act','hiv')]))
    hgbpcs.top = as.data.frame(table(hgbpcs[hgbpcs$pcs > hgbpcs.thr,c('act','hiv')]))
    hgbpcs.top$percent = hgbpcs.top$Freq/hgbpcs.tot$Freq*100

    figure.fn <- paste(fig.dir,'pcs_gene_bins_by_expr_and_hiv_barplot.pdf',sep='')
    cat('[8.1.2] Clustering of bins in act/silent genes with/out HIV -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(hgbpcs.top, aes(y=percent,x=act, fill=hiv)) +
        geom_bar(stat='identity',position='dodge',colour='black',alpha=0.6) +
        ylab(paste('Interactions with IHCD >=',hgbpcs.thr,'Hi-C reads/kb^2 [%]')) +
        xlab('Gene activity') +
        guides(fill=guide_legend(title='')) +
        theme_minimal()
    dev.off()

    ## The same now with SuperEnhancers

    hivbingen.df <- data.frame(chr = seqnames(hivbingen.gr),
                               beg = start(hivbingen.gr),
                               end = end(hivbingen.gr),
                               hiv = hivbingen.gr$hiv.cnt,
                               sen = hivbingen.gr$sen.close)

    hiv.gen.bin.file <- paste(hiv.gen.file,'sen.bins',sep='.')
    hiv.gen.bin.outf <- paste(hiv.gen.file,'sen.bins.out',sep='.')
    write.table(hivbingen.df,
                file      = hiv.gen.bin.file,
                quote     = FALSE,
                sep       ='\t',
                row.names = FALSE,
                col.names = TRUE
                )

    system(paste('python src/pairwise_interaction_bin_groups_SE.py',hiv.gen.bin.file, cool.file, hic.binsize, hiv.gen.bin.outf))

    hiv.gen.bin.pcs <- read.table(hiv.gen.bin.outf, sep='\t', col.names=c('chr.a','sen.a','hiv.a','span.a','chr.b','sen.b','hiv.b','span.b','pcs','hic.cnt'))

    hgbpcs <- hiv.gen.bin.pcs

    hgbpcs$sen <- NA
    hgbpcs$sen[hgbpcs$sen.a == hgbpcs$sen.b & hgbpcs$sen.a == 'SE'] = 'SE-SE'
    hgbpcs$sen[hgbpcs$sen.a == hgbpcs$sen.b & hgbpcs$sen.a == 'No SE'] = 'No SE-No SE'
    hgbpcs$sen[hgbpcs$sen.a != hgbpcs$sen.b] = 'SE-No SE'

    hgbpcs$sen = factor(hgbpcs$sen,c('SE-SE','SE-No SE','No SE-No SE'))

    hgbpcs$hiv <- NA
    hgbpcs$hiv[hgbpcs$hiv.a == hgbpcs$hiv.b & hgbpcs$hiv.a == 'HIV'] = 'HIV-HIV'
    hgbpcs$hiv[hgbpcs$hiv.a == hgbpcs$hiv.b & hgbpcs$hiv.a == 'No HIV'] = 'No HIV-No HIV'
    hgbpcs$hiv[hgbpcs$hiv.a != hgbpcs$hiv.b] = 'HIV-No HIV'

    hgbpcs.rev <- hgbpcs[,c('chr.b','sen.a','hiv.a','span.a','chr.a','sen.b','hiv.b','span.b','pcs','hic.cnt','sen','hiv')]
    names(hgbpcs.rev) <- c('chr.a','sen.a','hiv.a','span.a','chr.b','sen.b','hiv.b','span.b','pcs','hic.cnt','sen','hiv')
    hgbpcs.dup <- rbind(hgbpcs, hgbpcs.rev)

    ## Boxplot
    figure.fn <- paste(fig.dir,'pcs_gene_bins_by_senh_and_hiv.pdf',sep='')
    cat('[8.1.2] Clustering of bins in SE/No SE genes with/out HIV -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(data=hgbpcs, aes(x=sen, y=pcs, fill=hiv)) +
        geom_boxplot(outlier.shape=NA) +
        ylim(c(0,4e-4)) +
        ##    scale_y_log10() +
        theme_minimal() +
        guides(fill=guide_legend(title='')) +
        xlab('Super enhancer in gene proximity') +
        ylab('Interchromosomal Hi-C contact density [Hi-C reads/kb^2]')
    dev.off()

    ## Barplot
    hgbpcs.thr <-3e-4
    hgbpcs.tot = as.data.frame(table(hgbpcs[,c('sen','hiv')]))
    hgbpcs.top = as.data.frame(table(hgbpcs[hgbpcs$pcs > hgbpcs.thr,c('sen','hiv')]))
    hgbpcs.top$percent = hgbpcs.top$Freq/hgbpcs.tot$Freq*100

    figure.fn <- paste(fig.dir,'pcs_gene_bins_by_senh_and_hiv_barplot.pdf',sep='')
    cat('[8.1.2] Clustering of bins in SE/No SE genes with/out HIV -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=8, height=8)
    ggplot(hgbpcs.top, aes(y=percent,x=sen, fill=hiv)) +
        geom_bar(stat='identity',position='dodge',colour='black',alpha=0.6) +
        ylab(paste('Interactions with IHCD >=',hgbpcs.thr,'Hi-C reads/kb^2 [%]')) +
        xlab('Gene proximity to super enhancer') +
        guides(fill=guide_legend(title='')) +
        theme_minimal()
    dev.off()
}

##
## 8.2. 3D pattern of integrations in HIV-dense genes
##
cat('[8.2. 3D pattern of integrations in HIV-dense genes (HD genes)]\n')

## Sum HIV counts per gene.
rig.gr      <- genes.gr[genes.gr@ranges@width > hiv.gen.len]
rig.gr      <- rig.gr[grep("chr[0-9XY][0-9]*", as.vector(rig.gr@seqnames), value=FALSE)]
hiv.rig.ov  <- findOverlaps_nw(hiv.gr, rig.gr, ignore.strand = TRUE)
hiv.rig.cnt <- aggregate(hiv.cnt ~ gene.id,
                         data = data.frame(hiv.cnt = hiv.rig.ov@from, gene.id = hiv.rig.ov@to),
                         FUN=length)

rig.gr$hiv.cnt                      <- 0
rig.gr$hiv.cnt[hiv.rig.cnt$gene.id] <- hiv.rig.cnt$hiv.cnt
rig.gr$hiv.dens.kb                  <- rig.gr$hiv.cnt/rig.gr@ranges@width*1e3

## Define intervals
rig.pctl.intv    <- c(0.995,0.99,0.95,0.9,(sum(rig.gr$hiv.cnt==0)+1)/sum(rig.gr$hiv.cnt>=0),0.0)
rig.pctl.names   <- c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Sparse HIV','No HIV')
rig.pctl.colors  <- hue_pal()(length(rig.pctl.names))
## Sort intervals
rig.pctl.ord     <- order(rig.pctl.intv)
rig.pctl.intv.s  <- rig.pctl.intv[rig.pctl.ord]
rig.pctl.names.s <- rig.pctl.names[rig.pctl.ord]


## Annotate percentile position of genes in HIV-density
rig.gr$hiv.dens.pctl <- ecdf(rig.gr$hiv.dens.kb)(rig.gr$hiv.dens.kb)
rig.gr$hiv.dens.intv <- factor(rig.pctl.names.s[findInterval(rig.gr$hiv.dens.pctl,rig.pctl.intv.s)], rig.pctl.names)


## Barplot with hdgene intervals.
figure.fn <- paste(fig.dir,'hiv_hdgenes_density_interv_endogenous_expression.pdf',sep='')
cat('[8.2] Expression of endogenous genes vs HIV density (intervals boxplot) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(data=data.frame(intv=rig.gr$hiv.dens.intv, expr=rig.gr$expr), aes(x=intv, y=expr, fill=intv)) +
    geom_boxplot(outlier.shape=NA) +
    xlab('Genes by HIV density') +
    ylab('Endogenous expression [tpm]') +
    ylim(c(0,300))
dev.off()

## Figure expression-cluster-interval

## Add SE proximity info to rig.gr
rig.se.ov <- findOverlaps_nw(resize(rig.gr, width(rig.gr)+5000, fix='end'), sen.gr)
rig.gr$SE.close <- 'No superenhancer in gene proximity'
rig.gr$SE.close[unique(rig.se.ov@from)] <- 'Superenhancer in gene proximity'

## Add cluster info to intervals at gene TSS.
rig.tss.gr       <- resize(rig.gr,1)
rig.tss.gr.ov    <- findOverlaps_nw(rig.tss.gr, clust.gr)
rig.tss.gr$clust <- NA
rig.tss.gr$clust[rig.tss.gr.ov@from] <- clust.gr$global.cluster[rig.tss.gr.ov@to]

hiv.rig.intv.expr  <- data.frame(intv     = rig.tss.gr$hiv.dens.intv,
                                 hiv.dens = rig.tss.gr$hiv.dens.kb,
                                 clust    = rig.tss.gr$clust,
                                 expr     = rig.tss.gr$expr)
hiv.rig.intv.expr$intv <- as.character(hiv.rig.intv.expr$intv)
gen.de <- hiv.rig.intv.expr
hiv.rig.intv.expr[hiv.rig.intv.expr$intv %in% c("Top 0.5%", "Top 0.5%-1%", "Top 1%-5%"), 'intv'] <- "Top 5%"
hiv.rig.intv.expr$intv <- factor(hiv.rig.intv.expr$intv, c('Top 5%','Top 5%-10%','Sparse HIV','No HIV'))
hrie <- hiv.rig.intv.expr

figure.fn <- paste(fig.dir,'hiv_hdgenes_clust_intv_endogenous_expr.pdf',sep='')
cat('[8.2] Expression of endogenous genes vs HIV density and cluster -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(hrie[(hrie$intv %in% c('Sparse HIV','No HIV') | hrie$clust %in% c("A1","A2","B3")) & !is.na(hrie$clust),] ,aes(x=intv,y=expr,fill=clust)) +
    scale_fill_manual(name='3D Pattern at TSS', values = clust.colors) +
    geom_boxplot(outlier.shape=NA) +
    xlab('Genes by HIV density') +
    ylab('Endogenous gene expression [tpm]') +
    scale_y_log10(limits=c(1e-4, 1e3)) +
    theme_minimal()
dev.off()

## Stacked barplot, HIV dens vs gene expression.
agen.de <- gen.de[gen.de$expr > 0,]
sgen.de <- gen.de[gen.de$expr == 0,]
agen.de$expr.pctl <- ecdf(agen.de$expr)(agen.de$expr)
agen.de$expr.pctl <- as.factor(ceiling(agen.de$expr.pctl*10)/10)
sgen.de$expr.pctl <- as.factor('Silent')
gen.de            <- rbind(agen.de, sgen.de)
hgen.de           <- gen.de[gen.de$hiv.dens > 0,]
ngen.de           <- gen.de[gen.de$hiv.dens == 0,]
hgen.de$hiv.pctl  <- ecdf(hgen.de$hiv.dens)(hgen.de$hiv.dens)
hgen.de$hiv.pctl  <- as.factor(ceiling(hgen.de$hiv.pctl*10)/10)
ngen.de$hiv.pctl  <- as.factor('No HIV')
gen.de            <- rbind(hgen.de, ngen.de)
gen.de.cnt        <- as.data.frame(table(gen.de[,c('hiv.pctl','expr.pctl')]))
gen.de.cnt$expr.pctl <- factor(gen.de.cnt$expr.pctl, c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,'Silent'))
gen.de.cnt$hiv.pctl  <- factor(gen.de.cnt$hiv.pctl, c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,'No HIV'))
                                        #gen.de.cnt$intv      <- factor(gen.de.cnt$hiv.pctl, c("Top 0.5%", "Top 0.5%-1%", "Top 1%-5%", "Top 5%-10%", "Sparse HIV", "No HIV"))

## Plot HIV density vs gene Expression in 10 quantiles
figure.fn <- paste(fig.dir,'endogenous_expression_and_hiv_density.pdf',sep='')
cat('[8.2] Distribution of HIV density in endogenous expression quantiles -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=15, height=15)
ggplot(data=gen.de.cnt, aes(x=expr.pctl, y=Freq, fill=hiv.pctl)) +
    geom_bar(stat='identity', position='fill') +
    xlab('Endogenous gene expression (Percentile)') +
    ylab('Rate of HIV density in percentile groups') +
    guides(fill=guide_legend(title='HIV density (percentiles)')) +
    theme_minimal()
dev.off()

genes.hiv.clust.se <- data.frame(gene.name = rig.tss.gr$name,
                                 chr = seqnames(rig.tss.gr),
                                 beg = start(rig.gr),
                                 end = end(rig.gr),
                                 hiv.dens.kb = rig.tss.gr$hiv.dens.kb,
                                 expr.tpm = rig.tss.gr$expr,
                                 cluster = rig.tss.gr$clust,
                                 SE.close.5kb = as.factor(rig.tss.gr$SE.close)
                                 )

expr.w.hiv <- data.frame(hiv.dens = rig.tss.gr$hiv.dens.kb,
                         expr = rig.tss.gr$expr,
                         clust = rig.tss.gr$clust,
                         SE.close = as.factor(rig.tss.gr$SE.close)
                         )
                                        #expr.w.hiv <- expr.w.hiv[expr.w.hiv$expr > 1e-3 & expr.w.hiv$hiv.dens > 0,]
                                        #expr.w.hiv <- expr.w.hiv[expr.w.hiv$expr < 5e3 & expr.w.hiv$hiv.dens < 2,]
expr.w.hiv[expr.w.hiv$expr == 0, 'expr'] <- min(expr.w.hiv$expr)
expr.w.hiv[expr.w.hiv$hiv.dens == 0, 'hiv.dens'] <- min(expr.w.hiv$hiv.dens)
expr.w.hiv <- expr.w.hiv[!is.na(expr.w.hiv$clust),]
expr.w.hiv$clust <- factor(expr.w.hiv$clust,c('A1','A2','B1','B2','B3'))
## Extra plot: HIV dens vs Gene expression (by clusters):
figure.fn <- paste(fig.dir,'hiv_hdgenes_endogenous_expression_by_3dpatterns.pdf',sep='')
cat('[8.2] Expression of endogenous genes vs HIV density (by clusters) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=15, height=15)
ggplot(data=expr.w.hiv, aes(y=hiv.dens, x=expr, color=clust)) +
    geom_point(size=0.5) +
    geom_smooth(method=lm) +
    ##    scale_color_manual(name='', values=clust.colors) +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal() +
    ylab('HIV density [integ/kbp]') +
    xlab('Gene expression [tpm]')
dev.off()

## Extra plot: HIV dens vs Gene expression (by proximity to SE):
figure.fn <- paste(fig.dir,'hiv_hdgenes_endogenous_expression_by_proximity_to_SE.pdf',sep='')
cat('[8.2] Expression of endogenous genes vs HIV density (by SE proximity) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=15, height=15)
ggplot(data=expr.w.hiv, aes(y=hiv.dens, x=expr, color=SE.close)) +
    geom_point(size=0.5) +
    geom_smooth(method=lm) +
    ##    scale_color_manual(name='', values=clust.colors) +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal() +
    ylab('HIV density [integ/kbp]') +
    xlab('Gene expression [tpm]')
dev.off()

## Extra plot: Pie charts with genes close to SE, as function of HIV-target and Expr level
ewh = expr.w.hiv
ewh$hiv = 'No HIV'
ewh[ewh$hiv.dens > 0,]$hiv = 'With HIV'
ewh$exp.cat = 'Silent'
ewh.e = ewh[ewh$expr > 0,]
ewh.itv.names = c('Low 10%','Mid','High 10%')
ewh.e$exp.cat = factor(ewh.itv.names[findInterval(ecdf(ewh.e$expr)(ewh.e$expr), c(0,0.1,0.9))], ewh.itv.names)
ewh.tot = rbind(ewh[ewh$expr == 0,], ewh.e)
ewh.tot$exp.cat = factor(ewh.tot$exp.cat, c('Silent','Low 10%','Mid','High 10%'))

ggplot(data=as.data.frame(table(SE.close=ewh.tot$SE.close, hiv=ewh.tot$hiv, exp.cat=ewh.tot$exp.cat)), aes(x=hiv, y=Freq, fill=SE.close)) +
    geom_bar(stat='identity', position='fill') +
    xlab('') +
    ylab('') +
    coord_polar('y') +
    facet_wrap(~exp.cat)


## Add info to HD genes: distance to closest superEnhancer
rig.tss.gr@ranges <- rig.gr@ranges
hdg.sen.d <- distanceToNearest(rig.tss.gr, sen.gr, ignore.strand = TRUE)
rig.tss.gr$sen.d <- NA
rig.tss.gr$sen.d[hdg.sen.d@from] <- hdg.sen.d@elementMetadata$distance
rig.sen <- data.frame(clust = rig.tss.gr$clust,
                      intv  = rig.tss.gr$hiv.dens.intv,
                      sen.d = rig.tss.gr$sen.d)
rig.sen$intv = as.character(rig.sen$intv)
rig.sen[rig.sen$intv %in% c("Top 0.5%", "Top 0.5%-1%"), 'intv'] = "Top 1%"
rig.sen$intv = factor(rig.sen$intv, c('Top 1%','Top 1%-5%','Top 5%-10%','Sparse HIV','No HIV'))

## Figure: Distance to superEnhancers (grouped by intervals)
figure.fn <- paste(fig.dir,'hiv_hdgenes_intv_superenh.pdf',sep='')
cat('[8.2] Distance of HIV-dense to SuperEnhancers (grouped by intervals) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(data=rig.sen, aes(x=intv, y=sen.d, fill=intv)) +
    geom_boxplot(outlier.shape=NA) +
    scale_y_log10(limits=c(1e4,1e8)) +
    xlab('HIV-dense genes') +
    ylab('Distance to closest Super Enhancer [bp]')
dev.off()

rig.sen$intv = as.character(rig.sen$intv)
rig.sen[rig.sen$intv %in% c("Top 1%", "Top 1%-5%"), 'intv'] = "Top 5%"
rig.sen$intv = factor(rig.sen$intv, c('Top 5%','Top 5%-10%','Sparse HIV','No HIV'))


## Figure: Distance to superEnhancers (grouped by intervals and clusters)
figure.fn <- paste(fig.dir,'hiv_hdgenes_intv_clust_superenh.pdf',sep='')
cat('[8.2] Distance of HIV-dense to SuperEnhancers (grouped by intervals) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(data=rig.sen, aes(x=intv, y=sen.d, fill=clust)) +
    scale_fill_manual(name='3D pattern at TSS', values = clust.colors) +
    geom_boxplot(outlier.shape=NA) +
    scale_y_log10(limits=c(1e4,1e8)) +
    xlab('HIV-dense genes') +
    ylab('Distance to closest Super Enhancer [bp]')
dev.off()





## Overlap with HIV data
hiv.rig.intv.ov    <- findOverlaps_nw(hiv.gr, rig.gr, ignore.strand = TRUE)
hiv.rig.intv.table <- table(data.frame(intv  = rig.gr$hiv.dens.intv[hiv.rig.intv.ov@to],
                                       clust = factor(hiv.gr$clust[hiv.rig.intv.ov@from],
                                                      rev(c("A1","A2","B1","B2","B3")))
                                       )
                            )

hiv.rig.intv.rates <- as.data.frame(hiv.rig.intv.table/rowSums(hiv.rig.intv.table))

## Structure type: Stacked bar plot
figure.fn <- paste(fig.dir,'hiv_hdgenes_3dpatterns.pdf',sep='')
cat('[8.2] Distribution of HIV-dense genes in 3D patterns -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(hiv.rig.intv.rates, aes(x=intv, y=Freq, fill=clust)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name='', values = clust.colors) +
    xlab(paste('HIV-dense genes', sep='')) +
    ylab('3D pattern ratio of HIV integrations') +
    theme_minimal()
dev.off()

### Expressed/Silent genes: difference between targeted and non-targeted.
rig.tss.gr$expr.intv <- NA
rig.etss.gr <- rig.tss.gr[rig.tss.gr$expr > 0]
rig.etss.gr$expr.intv[rig.etss.gr$expr >= quantile(rig.etss.gr$expr,0.9)] <- 'Top 10%'
rig.etss.gr$expr.intv[rig.etss.gr$expr <= quantile(rig.etss.gr$expr,0.1)] <- 'Bottom 10%'
rig.tss.gr$expr.intv[rig.tss.gr$expr == 0] <- 'Silent'
rig.tss.expr.gr <- c(rig.tss.gr[rig.tss.gr$expr == 0], rig.etss.gr[rig.etss.gr$expr.intv %in% c('Top 10%','Bottom 10%')])
rig.tss.expr.gr$expr.intv <- as.factor(rig.tss.expr.gr$expr.intv)
rig.tss.expr.gr$hiv.target <- NA
rig.tss.expr.gr$hiv.target[rig.tss.expr.gr$hiv.cnt > 0] <- 'HIV Target'
rig.tss.expr.gr$hiv.target[rig.tss.expr.gr$hiv.cnt == 0] <- 'No HIV'
rig.tss.expr.gr$hiv.target <- as.factor(rig.tss.expr.gr$hiv.target)

rig.tss.df <- data.frame(expr.intv  = rig.tss.expr.gr$expr.intv,
                         clust      = rig.tss.expr.gr$clust,
                         senh.dist  = rig.tss.expr.gr$sen.d,
                         hiv.target = rig.tss.expr.gr$hiv.target,
                         hiv.dens   = rig.tss.expr.gr$hiv.dens.kb,
                         gene.expr  = rig.tss.expr.gr$expr
                         )

rig.tss.df$expr.intv = factor(rig.tss.df$expr.intv, c('Top 10%', 'Bottom 10%', 'Silent'))
## Pie chart of HIV targets in Silent genes and Top 10% expressed genes.
figure.fn <- paste(fig.dir,'expressed_genes_hiv_targets.pdf',sep='')
cat('[8.2] Percentage of HIV targets among expressed and silent genes -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=5)
ggplot(data=as.data.frame(table(rig.tss.df[,c('expr.intv','hiv.target')])), aes(x="", y=Freq, fill=hiv.target)) +
    geom_bar(width = 1, stat = "identity", position="fill") +
    coord_polar('y') +
    xlab('') +
    ylab('') +
    guides(fill=guide_legend(title='')) +
    facet_wrap(~expr.intv)
dev.off()

## Pie chart of 3D pattern in Silent genes and Top 10% expressed genes.
figure.fn <- paste(fig.dir,'expressed_genes_hiv_3dpattern.pdf',sep='')
cat('[8.2] 3D patterns of HIV targets and non-targets in expressed and silent genes -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(data=as.data.frame(table(rig.tss.df[,c('expr.intv','hiv.target','clust')])), aes(x="", y=Freq, fill=clust)) +
    geom_bar(width = 1, stat = "identity", position="fill") +
    scale_fill_manual(name='3D Pattern', values = clust.colors) +
    coord_polar('y') +
    xlab('') +
    ylab('') +
    facet_wrap(~expr.intv+hiv.target, ncol=2)
dev.off()



## Distance to SuperEnhancers.
figure.fn <- paste(fig.dir,'expressed_genes_hiv_dist_to_superenhancer.pdf',sep='')
cat('[8.2] Distance to closest SuperEnhancer of HIV targets and non-targets in expressed and silent genes -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=4, height=8)
ggplot(data=rig.tss.df, aes(x=expr.intv, y=senh.dist, fill=hiv.target)) +
    geom_boxplot(outlier.size=0.3) +
    scale_y_log10() +
    xlab('') +
    ylab('')
dev.off()

## Same with full gene bodies (No TSS and no cluster info)
## rig.gr <- rig.gr[rig.gr$expr > 0]
## rig.gr$expr.intv <- NA
## rig.gr$expr.intv[rig.gr$expr >= quantile(rig.gr$expr,0.9)] <- 'Top 10%'
## rig.gr$expr.intv[rig.gr$expr <= quantile(rig.gr$expr,0.1)] <- 'Bottom 10%'
## rig.expr.gr <- rig.gr[rig.gr$expr.intv %in% c('Top 10%','Bottom 10%')]
## rig.expr.gr$expr.intv <- factor(rig.expr.gr$expr.intv, c('Top 10%','Bottom 10%'))
## rig.expr.gr$hiv.target <- NA
## rig.expr.gr$hiv.target[rig.expr.gr$hiv.cnt > 0] <- 'HIV Target'
## rig.expr.gr$hiv.target[rig.expr.gr$hiv.cnt == 0] <- 'No HIV'
## rig.expr.gr$hiv.target <- as.factor(rig.expr.gr$hiv.target)
## Compute contacts between SuperEnhancer and Gene TSS
rig.expr.gr <- rig.tss.expr.gr

rig.expr.infile <- paste(hiv.gen.file,'hiv.targets',sep='.')
rig.expr.outf   <- paste(hiv.gen.file,'hiv.targets.out',sep='.')
rig.expr.ipcs   <- data.frame(name.1=c(), dens=c(), act=c(), sen.id=c(), sen.1=c(), sen.2=c(), pcs=c(), reads=c(), gen.exp=c(), hiv.target=c())
rig.expr.pcs    <- data.frame(name.1=c(), dens=c(), act=c(), sen.id=c(), sen.1=c(), sen.2=c(), pcs=c(), reads=c(), gen.exp=c(), hiv.target=c())
for (gen.exp in unique(rig.expr.gr$expr.intv)) {
    for (hiv.target in unique(rig.expr.gr$hiv.target)) {
        rig.tmp.gr <- rig.expr.gr[rig.expr.gr$expr.intv == gen.exp & rig.expr.gr$hiv.target == hiv.target]
        tmp.df <- data.frame(name    = rig.tmp.gr$name,
                             chr     = seqnames(rig.tmp.gr),
                             beg     = rig.tmp.gr@ranges@start,
                             end     = rig.tmp.gr@ranges@start + rig.tmp.gr@ranges@width,
                             dens.kb = rig.tmp.gr$hiv.dens.kb,
                             act     = rig.tmp.gr$active
                             )
        ## Write gene file for processing.
        write.table(tmp.df,
                    file      = rig.expr.infile,
                    quote     = FALSE,
                    sep       ='\t',
                    row.names = FALSE,
                    col.names = TRUE
                    )
        
        ## Run Group PCS (intrachromosomal)
        cat('[8.2] Computing intrachromosomal contact scores between SuperEnhancers and',gen.exp,'with condition:',hiv.target,'\n')
        system(paste('python src/pairwise_intra_interaction_score_group.py',rig.expr.infile, hp.sen.fname, cool.file, hic.binsize, rig.expr.outf))

        out.df <- read.table(rig.expr.outf, col.names = c('name.1','dens','act','sen.id','sen.1','sen.2','pcs','reads'))
        out.df$expr.intv  <- gen.exp
        out.df$hiv.target <- hiv.target
        rig.expr.pcs      <- rbind(rig.expr.pcs, out.df)
        
        ## Run Group PCS (interchromosomal)
        cat('[8.2] Computing interchromosomal contact scores between SuperEnhancers and',gen.exp,'with condition:',hiv.target,'\n')
        ##        system(paste('python src/pairwise_interaction_score_group.py',rig.expr.infile, hp.sen.fname, cool.file, hic.binsize, rig.expr.outf))

        out.df <- read.table(rig.expr.outf, col.names = c('name.1','dens','act','sen.id','sen.1','sen.2','pcs','reads'))
        out.df$expr.intv  <- gen.exp
        out.df$hiv.target <- hiv.target
        rig.expr.ipcs     <- rbind(rig.expr.ipcs, out.df)
    }
}

## Store data
write.table(rig.expr.pcs,
            file = 'annotations/hiv-target-superenh-pcs.txt',
            quote     = FALSE,
            sep       ='\t',
            row.names = FALSE,
            col.names = TRUE
            )

write.table(rig.expr.ipcs,
            file = 'annotations/hiv-target-superenh-ipcs.txt',
            quote     = FALSE,
            sep       ='\t',
            row.names = FALSE,
            col.names = TRUE
            )


## Density plot
figure.fn <- paste(fig.dir,'expressed_genes_tss_senh_intrachr_pcs.pdf',sep='')
cat('[9.1] Contact score between HD genes and superenhancers -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(rig.expr.max.pcs, aes(y=..density..,x=pcs, colour=expr.intv, linetype=hiv.target)) +
    geom_freqpoly(binwidth=0.005) +
    scale_y_log10() +
    xlim(c(0,0.1)) +
    xlab('Intrachromosomal pairwise contact score [Hi-C reads/kb^2]') +
    ylab('Density') +
    guides(colour=guide_legend(title='')) +
    ggtitle('Clustering with SuperEnhancers (Jurkat)') +
    theme_minimal()
dev.off()

## ecdf
figure.fn <- paste(fig.dir,'expressed_genes_tss_senh_intrachr_pcs_ecdf.pdf',sep='')
cat('[9.1] Contact score between HD genes and superenhancers (ecdf) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(rig.expr.pcs, aes(x=pcs, linetype=hiv.target, colour=expr.intv)) +
    ##    stat_ecdf(geom='step') +
    stat_ecdf() +
    ##    scale_y_log10() +
    xlim(c(0,0.05)) +
    xlab('Intrachromosomal pairwise contact score [Hi-C reads/kb^2]') +
    ylab('ecdf') +
    guides(colour=guide_legend(title='')) +
    ggtitle('Clustering with SuperEnhancers (Jurkat)') +
    theme_minimal()
dev.off()

## Filter only the maximum PCS value for each gene.
rig.expr.max.pcs <- aggregate(pcs~name.1+expr.intv+hiv.target, data=rig.expr.pcs, FUN=max)

## PCS ecdf
figure.fn <- paste(fig.dir,'expressed_genes_tss_senh_intrachr_max_pcs_ecdf.pdf',sep='')
cat('[9.1] Contact score between HD genes and superenhancers (ecdf) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(rig.expr.max.pcs, aes(x=pcs, linetype=hiv.target, colour=expr.intv)) +
    ##    stat_ecdf(geom='step') +
    stat_ecdf() +
    ##    scale_y_log10() +
    xlim(c(0,1)) +
    xlab('Intrachromosomal pairwise contact score [Hi-C reads/kb^2]') +
    ylab('ecdf') +
    guides(colour=guide_legend(title=''), linetype=guide_legend(title='')) +
    ggtitle('Gene contacts with SuperEnhancers (Jurkat)') +
    theme_minimal()
dev.off()

## Filter only the maximum PCS value for each gene.
rig.expr.max.ipcs <- aggregate(pcs~name.1+expr.intv+hiv.target, data=rig.expr.ipcs, FUN=max)

## IPCS ecdf
figure.fn <- paste(fig.dir,'expressed_genes_tss_senh_intrachr_max_ipcs_ecdf.pdf',sep='')
cat('[9.1] Contact score between HD genes and superenhancers (ecdf) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(rig.expr.max.ipcs, aes(x=pcs, linetype=hiv.target, colour=expr.intv)) +
    ##    stat_ecdf(geom='step') +
    stat_ecdf() +
    ##    scale_y_log10() +
    xlim(c(0,0.01)) +
    xlab('Interchromosomal pairwise contact score [Hi-C reads/kb^2]') +
    ylab('ecdf') +
    guides(colour=guide_legend(title='')) +
    ggtitle('Clustering with SuperEnhancers (Jurkat)') +
    theme_minimal()
dev.off()


## Barplot
rig.sen.pcs.thr  <- 0.005
rig.sen.pcs.data <- as.data.frame(table(rig.expr.pcs[rig.expr.pcs$pcs > rig.sen.pcs.thr, c('expr.intv','hiv.target')]))
rig.sen.pcs.data.all <- as.data.frame(table(rig.expr.pcs[, c('expr.intv','hiv.target')]))
rig.sen.pcs.data$Freq <- rig.sen.pcs.data$Freq/rig.sen.pcs.data.all$Freq
ggplot(rig.sen.pcs.data, aes(y=Freq,x=expr.intv, fill=hiv.target)) +
    geom_bar(stat='identity',colour='black',alpha=0.6,position='dodge') +
    ylab(paste('Interactions with IPCS >=',rig.sen.pcs.thr,'reads/kb^2 [%]')) +
    xlab('') +
    ggtitle('Contacts with SuperEnhancers') +
    theme_minimal()

##
## 8.3. Chromatin signal in HD gene bodies.
##
cat('[8.3. Chromatin signal in HD gene TSS and TES]\n')
cm.raw.path = 'chip_datasets/raw_select/*.bed'

## Load data from raw Chip files
cat('[8.3] Read raw ChIP-seq data from:',cm.raw.path,'\n')
cm.raw.gr    = GRanges()
cm.raw.names = c()

## Minimum gene length for analysis.
rig.gen.len  = rig.gen.tss + rig.gen.tes

## Convert Top 0.5% to 1%
hdg.gr = rig.gr
hdg.gr = hdg.gr[hdg.gr@ranges@width >= rig.gen.len]
hdg.gr$hiv.dens.intv = as.character(hdg.gr$hiv.dens.intv)
hdg.gr[hdg.gr$hiv.dens.intv == 'Top 0.5%']$hiv.dens.intv = 'Top 1%'
hdg.gr[hdg.gr$hiv.dens.intv == 'Top 0.5%-1%']$hiv.dens.intv = 'Top 1%'

dataset.id = 1
for (f in Sys.glob(cm.raw.path)) {
    cat('[8.3] reading file: ',f,'...\n',sep='')
    ## Extract CM name from filename (text before '_')
    path         <- strsplit(f,'/')[[1]]
    fname        <- strsplit(path[length(path)],'.',fixed=T)[[1]][1]
    cm.name      <- strsplit(fname,'_',fixed=T)[[1]][1]
    cm.raw.names <- c(cm.raw.names,cm.name)
    ## Read file (only chr[1-9XY][0-9]*)
    d            <- subset(read.table(f), grepl("^chr[0-9XY][0-9]*$",V1))
    d$V1         <- factor(d$V1, unique(d$V1))
    cm.tmp.gr    <- GRanges(seqnames = Rle(d$V1),
                            ranges   = IRanges(start=d$V2,end=d$V3),
                            score    = d$V4,
                            dataset  = rep(dataset.id,nrow(d)),
                            name     = rep(cm.name,nrow(d)),
                            )
    ## Fill spaces with score 0
    cm.pad.gr         <- setdiff(genome.gr, cm.tmp.gr)
    cm.pad.gr$score   <- 0
    cm.pad.gr$dataset <- dataset.id
    cm.pad.gr$name    <- cm.name
    ## Append to global GRanges
    cm.raw.gr    <- c(cm.raw.gr, cm.tmp.gr, cm.pad.gr)
    ## Verbose dataset info
    cat('[8.3] raw chip-seq #',dataset.id,': ',cm.name,' (',nrow(d),' lines)\n',sep='')
    ## Increase loop counter
    dataset.id <- dataset.id + 1
}

hdg.signal.data = data.frame(x=c(), clu=c(), intv=c(), chip=c(), anchor=c(), score=c())
for (intv in unique(hdg.gr$hiv.dens.intv)) {
    cat('[8.3] generating GRanges intervals for ',intv,'...\n',sep='')
    tmp.gr.p <- hdg.gr[hdg.gr$hiv.dens.intv == intv & hdg.gr@strand == '+']
    tmp.gr.n <- hdg.gr[hdg.gr$hiv.dens.intv == intv & hdg.gr@strand == '-']
    intv.gr <- GRanges()
    ## TSS
    for (x in seq(-rig.gen.tss, rig.gen.tss, rig.gen.step)) {
        ## Create interval around gene TSS
        intv.gr <- c(intv.gr,
                     GRanges(seqnames = Rle(tmp.gr.p@seqnames),
                             ranges   = IRanges(start=tmp.gr.p@ranges@start + x, end=tmp.gr.p@ranges@start + x + rig.gen.bin),
                             x        = rep(x, length(tmp.gr.p@seqnames)),
                             anchor   = rep('TSS', length(tmp.gr.p@seqnames)),
                             active   = tmp.gr.p$active
                             ),
                     GRanges(seqnames = Rle(tmp.gr.n@seqnames),
                             ranges   = IRanges(start=tmp.gr.n@ranges@start + tmp.gr.n@ranges@width - x - rig.gen.bin, end=tmp.gr.n@ranges@start + tmp.gr.n@ranges@width - x),
                             x        = rep(x, length(tmp.gr.n@seqnames)),
                             anchor   = rep('TSS', length(tmp.gr.n@seqnames)),
                             active   = tmp.gr.n$active
                             )
                     )
    }

    ## TES
    for (x in seq(-rig.gen.tes, rig.gen.tes, rig.gen.step)) {
        ## Create interval around gene TES
        intv.gr <- c(intv.gr,
                     GRanges(seqnames = Rle(tmp.gr.p@seqnames),
                             ranges   = IRanges(start=tmp.gr.p@ranges@start + tmp.gr.p@ranges@width + x, end=tmp.gr.p@ranges@start + tmp.gr.p@ranges@width + x + rig.gen.bin),
                             x        = rep(x, length(tmp.gr.p@seqnames)),
                             anchor   = rep('TES', length(tmp.gr.p@seqnames)),
                             active   = tmp.gr.p$active
                             ),
                     GRanges(seqnames = Rle(tmp.gr.n@seqnames),
                             ranges   = IRanges(start=tmp.gr.n@ranges@start - x - rig.gen.bin, end=tmp.gr.n@ranges@start - x),
                             x        = rep(x, length(tmp.gr.n@seqnames)),
                             anchor   = rep('TES', length(tmp.gr.n@seqnames)),
                             active   = tmp.gr.n$active
                             )
                     )
    }
    ## Add 3D pattern information to TSS/TES intervals
    intv.clu.ov  <- findOverlaps_nw(intv.gr, clust.gr, ignore.strand = TRUE)
    intv.clu.ref <- data.frame(intv.id = intv.clu.ov@from, clu = clust.gr[intv.clu.ov@to]$global.cluster)
    intv.gr$clust <- NA
    intv.gr@elementMetadata$clust[intv.clu.ref$intv.id] <- as.character(intv.clu.ref$clu)

    ## Overlap with data
    tmp.ov  <- findOverlaps_nw(cm.raw.gr, intv.gr, ignore.strand = TRUE)
    
    hdg.signal.data <- rbind(hdg.signal.data,
                             data.frame(x      = intv.gr$x[tmp.ov@to],
                                        clu    = intv.gr$clust[tmp.ov@to],
                                        intv   = intv,
                                        chip   = cm.raw.gr$name[tmp.ov@from],
                                        anchor = intv.gr$anchor[tmp.ov@to],
                                        score  = cm.raw.gr$score[tmp.ov@from],
                                        active = intv.gr$active[tmp.ov@to]
                                        )
                             )
}

## Aggregate signal to compute median and 5-95 percentiles.
hdg.signal.agg = do.call(data.frame, aggregate(score ~ x+intv+chip+anchor+active, data = hdg.signal.data, FUN = function(x) c(median = median(x), mean = mean(x))))

hdg.signal.clu.agg = do.call(data.frame, aggregate(score ~ x+chip+clu+anchor+active, data = hdg.signal.data, FUN = function(x) c(median = median(x), mean = mean(x))))

hdg.signal.all = do.call(data.frame, aggregate(score ~ x+chip+intv+clu+anchor+active, data = hdg.signal.data, FUN = function(x) c(median = median(x), mean = mean(x))))

### HOTSPOTS
## Sort TSS/TES
hdg.signal.agg$anchor = factor(hdg.signal.agg$anchor, c('TSS','TES'))
hdg.signal.agg$active = as.character(hdg.signal.agg$active)
hdg.signal.agg[hdg.signal.agg$active == TRUE,'active'] = 'Active'
hdg.signal.agg[hdg.signal.agg$active == FALSE,'active'] = 'Silent'
hdg.signal.agg$active = factor(hdg.signal.agg$active, c('Active','Silent'))

## Figure with TSS/TES signal for each chip-seq profile
figure.fn <- paste(fig.dir,'hiv_hdgenes_tss_tes_chip_hotspots.pdf',sep='')
cat('[8.3] Chromatin signal in HD gene TSS and TES plots (hotspots) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=14, height=14)
ggplot(data=hdg.signal.agg[hdg.signal.agg$intv == c('No HIV','Top 1%','Top 5%-10%'),], aes(x=x, y=score.mean, color=intv, linetype=active)) +
    geom_line() +
    facet_wrap(~chip+anchor, scales='free_y')
dev.off()

## 3D PATTERNS
## Sort TSS/TES
hdg.signal.clu.agg$anchor = factor(hdg.signal.clu.agg$anchor, c('TSS','TES'))
hdg.signal.clu.agg$active = as.character(hdg.signal.clu.agg$active)
hdg.signal.clu.agg[hdg.signal.clu.agg$active == TRUE,'active'] = 'Active'
hdg.signal.clu.agg[hdg.signal.clu.agg$active == FALSE,'active'] = 'Silent'
hdg.signal.clu.agg$active = factor(hdg.signal.clu.agg$active, c('Active','Silent'))

## Figure with TSS/TES signal for each chip-seq profile
figure.fn <- paste(fig.dir,'hiv_hdgenes_tss_tes_chip_clust.pdf',sep='')
cat('[8.3] Chromatin signal in HD gene TSS and TES plots (clusters) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=14, height=14)
ggplot(data=hdg.signal.clu.agg, aes(x=x, y=score.mean, color=clu, linetype=active)) +
    geom_line() +
    facet_wrap(~chip+anchor, scales='free_y')
dev.off()

## BOTH
hdg.signal.all$anchor = factor(hdg.signal.all$anchor, c('TSS','TES'))
hdg.signal.all$active = as.character(hdg.signal.all$active)
hdg.signal.all[hdg.signal.all$active == TRUE,'active'] = 'Active'
hdg.signal.all[hdg.signal.all$active == FALSE,'active'] = 'Silent'
hdg.signal.all$active = factor(hdg.signal.all$active, c('Active','Silent'))

## Figure with TSS/TES signal for each chip-seq profile
figure.fn <- paste(fig.dir,'hiv_hdgenes_tss_tes_chip_all.pdf',sep='')
cat('[8.3] Chromatin signal in HD gene TSS and TES plots (hotspots and clusters)-> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=14, height=14)
ggplot(data=hdg.signal.all[hdg.signal.all$intv == c('No HIV','Top 1%','Top 5%-10%'),], aes(x=x, y=score.mean, color=clu, linetype=active)) +
    geom_line() +
    facet_wrap(~chip+intv+anchor, scales='free_y')
dev.off()

## Genome browser-like figure.

## Read clusters in bed format, reassign cluster order. (Update with each cluster computation)
##clust.bed = read.table(cluster.list[[def.clust]], header=F, col.names=c('chrom','start','end','name'))
clust.bed = data.frame(chrom=seqnames(clust.gr), start=start(clust.gr), end=end(clust.gr), name=clust.gr$global.cluster)
clust.bed$score = 1 #B3
clust.bed$score[clust.bed$name == 'A1'] = 5 #A1
clust.bed$score[clust.bed$name == 'B1'] = 3 #B1
clust.bed$score[clust.bed$name == 'B2'] = 2 #B2
clust.bed$score[clust.bed$name == 'A2'] = 4 #A2

## HIV, genes and SE datasets in BED format.
hiv.bed = data.frame(chrom=hiv.gr@seqnames, start=hiv.gr@ranges@start, end=hiv.gr@ranges@start+1)
sen.bed = data.frame(chrom=sen.data$chr, start=sen.data$beg, end=sen.data$end)
gen.gr = sort(genes.gr)
gen.bed = data.frame(chrom=seqnames(gen.gr), start=start(gen.gr), end=end(gen.gr), value=gen.gr$expr, act=gen.gr$active)
gen.bed$row = 1
gen.bed$row[gen.bed$act] = 2

## Define Figure regions.
chr  = c('chr16','chr2','chr1','chr1','chr1','chr17','chr19','chr5')
cbeg = c(73e6,96e6,40e6,157e6,225e6,62e6,5e6,50e6)
cend = c(77e6,100e6,44e6,161e6,235e6,66e6,30e6,60e6)
hcol = c(100,100,100,100,200,100,300,200)

signals = c('H3K27Ac','H3K9me3','CTCF','BRD4','MYC')

## Iterate over all Figure regions.
for (i in seq(1,length(chr))) {
    ## Call Python script, store to file (autogen name).
    hic_coords = paste(chr[i],':',formatC(cbeg[i],format='d'),'-',formatC(cend[i],format='d'),sep='')
    cat('Fetching matrix from Hi-C file for region:',hic_coords,'\n')
    system(paste('python src/cooler_to_mat.py', cool.file, hic_coords,'annotations/hic_tmp_mat.txt'))
    ## Read Hi-C matrix.
    hic_m = as.matrix(read.table('annotations/hic_tmp_mat.txt', header=T, row.names=1, check.names=F))

    ## Figure with TSS/TES signal for each chip-seq profile (autogen name based on region).
    figure.fn <- paste(fig.dir,'genome_',chr[i],'_',formatC(cbeg[i],format='d'),'_',formatC(cend[i],format='d'),'.pdf',sep='')
    cat('[8.3] Genome browser-like figure ',chr[i],':',cbeg[i],'-',cend[i],' -> ',figure.fn,'\n',sep='')
    pdf(figure.fn, useDingbats=F, width=14, height=14)
    sigcnt = length(signals)
    layout(matrix(seq(1,sigcnt+5),sigcnt+5,1,byrow=T), widths=c(1), heights=c(c(10,3,5,1,2),rep(2,sigcnt-1),4))
    par(mar=c(1,3.5,1,2))
    phic = plotHic(log1p(hic_m), chr[i], cbeg[i], cend[i], max_y = hcol[i], zrange=c(0,6), palette = SushiColors(7))
    par(mar=c(0,3.5,0,2))
    plotBed(beddata = hiv.bed, chrom=chr[i], chromstart=cbeg[i], chromend=cend[i], ylab='HIV')
    mtext('HIV',side=2,line=1.75,cex=0.6,font=2)
                                        #mtext('HIV',side=2,line=1.75,cex=1,font=2)
    par(mar=c(0,3.5,0,2))
    clust.bed.allrows = rbind(clust.bed,data.frame(chrom=rep(chr[i],5),start=rep(cbeg[i],5),end=rep(cbeg[i]+1),name=rep('B3',5),score=c(1,2,3,4,5)))
    plotBed(beddata = clust.bed.allrows, chrom=chr[i], chromstart=cbeg[i], chromend=cend[i], colorby=clust.bed.allrows$score, colorbycol=SushiColors(5), row='given',rownumber=clust.bed.allrows$score,type='region',rowlabels=rev(c('A1','A2','B1','B2','B3')),rowlabelcex=0.75, plotbg='grey95')
    par(mar=c(0,3.5,0,2))
    plotBed(beddata = sen.bed, chrom=chr[i], chromstart=cbeg[i], chromend=cend[i], rowlabels='SE', rowlabelcex=0.75, rownumber=1, row='given')
    par(mar=c(0,3.5,0,2))
    plotBed(beddata = gen.bed, chrom=chr[i], chromstart=cbeg[i], chromend=cend[i], rowlabelcex=0.75, row='given', rownumber=gen.bed$row, type='region', rowlabels=c('silgen', 'actgen'))
    for (signame in signals) {
        cr.sub = cm.raw.gr[cm.raw.gr@seqnames == chr[i] & cm.raw.gr@ranges@start >= cbeg[i] & cm.raw.gr@ranges@start < cend[i] & cm.raw.gr$name == signame]
        bg.data = data.frame(chrom=cr.sub@seqnames, start=cr.sub@ranges@start, end=cr.sub@ranges@start+cr.sub@ranges@width-1, value=cr.sub$score)
        if (signame == signals[length(signals)]) {
            par(mar=c(3.5,3.5,0,2))
            plotBedgraph(bg.data, chr[i], cbeg[i], cend[i])
            #axis(side=2,las=2,cex=0.5,font=2)
            mtext(signame,side=2,line=1.75,cex=0.6,font=2)
            labelgenome(chr[i],cbeg[i],cend[i],scale='Mb')
        } else {
            par(mar=c(0,3.5,0,2))
            plotBedgraph(bg.data, chr[i], cbeg[i], cend[i])
            #axis(side=2,las=2,cex=0.5,font=2)
            mtext(signame,side=2,line=1.75,cex=0.6,font=2)
        }
    }
    dev.off()
}





##
## 9. SuperEnhancers and HD genes.
##
cat('[9. SuperEnhancers and HD genes]\n')

##
## 9.1. Compute PCS between HD genes and SuperEnhancers (DBsuper)
##
cat('[9.1. PCS between HD genes and DBsuper]\n')

## Output file
sen.gen.outf  = paste(sen.pcs.file,'gen.out',sep='.')

## Run Group PCS
system(paste('python src/pairwise_interaction_score_group.py',hiv.gen.file, hp.sen.fname, cool.file, hic.binsize, sen.gen.outf))

sen.gen.pcs <- read.table(sen.gen.outf, col.names = c('name.1','dens','act','sen.id','sen','sen.id','pcs','reads'))

sen.gen.pcs$hiv = 'No HIV'
sen.gen.pcs[sen.gen.pcs$dens > 0,'hiv'] = 'RIG'

## Plot
figure.fn <- paste(fig.dir,'hiv_hdgenes_superenhancer.pdf',sep='')
cat('[9.1] Contact score between HD genes and superenhancers -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(sen.gen.pcs, aes(y=..density..,x=pcs, colour=hiv, linetype=act)) +
    geom_freqpoly(binwidth=0.001) +
    scale_y_log10() +
    xlim(c(0,0.01)) +
    xlab('Interchromosomal pairwise contact score [Hi-C reads/kb^2]') +
    ylab('Density') +
    guides(colour=guide_legend(title='')) +
    ggtitle('Clustering of HIV-dense genes with SuperEnhancers (Jurkat)') +
    theme_minimal()
dev.off()

##
## 9.2. Comparison of SuperEnhancers by Cluster
##
sen.extend.region = 5e3
sen.clu.ov = findOverlaps_nw(sen.gr + sen.extend.region, clust.gr, ignore.strand=TRUE)
sen.clu.idx = data.frame(sen.id = sen.clu.ov@from, clu = clust.gr[sen.clu.ov@to]$global.cluster)
sen.clu.idx = aggregate(clu ~ sen.id, data = sen.clu.idx, FUN=MaxTable)
sen.gr$clust = NA
sen.gr$clust[sen.clu.idx$sen.id] = sen.clu.idx$clu
sen.hiv.ov = findOverlaps_nw(sen.gr + sen.extend.region, hiv.gr, ignore.strand=TRUE)
sen.hiv.cnt = aggregate(cnt ~ sen.id, data= data.frame(sen.id = sen.hiv.ov@from, cnt=sen.hiv.ov@to), FUN=length)
sen.gr$hiv.cnt = 0
sen.gr$hiv.cnt[sen.hiv.cnt$sen.id] = sen.hiv.cnt$cnt
sen.gr$hiv.dens.kb = sen.gr$hiv.cnt / sen.gr@ranges@width * 1e3

hiv.sen.df = data.frame(clust = sen.gr$clust, width = sen.gr@ranges@width, hiv.cnt = sen.gr$hiv.cnt)
hiv.sen.agg = aggregate(. ~ clust, data=hiv.sen.df, FUN=sum)
hiv.sen.agg$hiv.dens.mb = hiv.sen.agg$hiv.cnt / hiv.sen.agg$width * 1e6

## Barplot of categories grouped by 3d pattern.
figure.fn <- paste(fig.dir,'hiv_density_superenhancer_clust.pdf',sep='')
cat('[4.5] HIV integration density (groups: SuperEnhancer) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=5,height=8)
ggplot(data=hiv.sen.agg, aes(fill=clust,x=clust,y=hiv.dens.mb)) +
    geom_bar(stat='identity',position='dodge',colour='black') +
    scale_fill_manual(name='Clusters',values=c(clust.colors)) +
    xlab('3D patterns') +
    ylab('HIV integrations in SuperEnhancers per Mbp')
dev.off()


