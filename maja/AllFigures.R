load("genidf.Robj")
library(ggplot2)
genidf$RIGSLISTE <- ifelse(genidf$numberOfLists>=2, "2 or more", 
                           ifelse(genidf$numberOfLists==1, "1", "0"))



SupplementaryFig1D <- ggplot(genidf, aes(x=RIGSLISTE, fill=SuperenhancersA))+
     geom_bar(position="fill")+
     ylab("")+
     xlab("")+
     theme_light()+
     coord_polar(theta="y")+
     scale_fill_discrete("RIGs lists")+
     theme(legend.position = "bottom")
ggsave(SupplementaryFig1D, file="SupplementaryFigure1D.pdf", width=5, height = 5)

Fig1C <- ggplot(genidf, aes(x=factor(numberOfLists), y=distanceToNearestSE))+
     stat_boxplot(geom="errorbar", coef=2.3 )+
     geom_boxplot(outlier.color = NA)+
     ylab("Distance to nearest SE (bp)")+
     xlab("Number of datasets with HIV-1 target gene")+
     theme_light()+
     scale_fill_discrete("RIGs lists")+coord_cartesian(ylim=c(0,5e06))
ggsave(Fig1C, file="Figure1C.pdf", width=7, height = 5)


gg <- data.frame(i.Mean=genidf[,c("i.Mean")],
                 Mean=genidf[,c("Mean")],
                 numberOfLists = factor(ifelse(genidf$numberOfLists>0, 
                                               "Genes with integrations",
                                               "Genes withouth integrations"), 
                                        levels=c("Genes withouth integrations", 
                                                 "Genes with integrations")),
                 gene_biotype=genidf$gene_biotype
                 
)
Fig2A <- ggplot(gg[gg$gene_biotype=="protein_coding",], aes(numberOfLists,i.Mean))+
     geom_violin()+
     scale_y_continuous()+
     ylab("Rlog of mean expression value")+
     theme_light()+
     xlab("Activated control cells")
ggsave(Fig2A, file="Figure2A.pdf", width=7, height = 5)
Fig2B <- ggplot(genidf[gene_biotype=="protein_coding"],aes(x=as.factor(numberOfLists),
                           group=numberOfLists,
                           y=i.Mean))+
     stat_boxplot(geom="errorbar")+
     geom_boxplot(outlier.colour = NA)+
     ylab("Rlog of mean expression in activated cells")+
     xlab("Number of lists")+
     theme_light()
ggsave(Fig2B, file="Figure2B.pdf", width=7, height = 5)
Fig2C <- ggplot(genidf[gene_biotype=="protein_coding"], aes(x=RIGSLISTE, fill=RIGSLISTE,
                                                            y= i.Mean))+
     stat_boxplot(geom="errorbar")+
     geom_boxplot(outlier.color = NA)+
     facet_wrap(~SuperenhancersA)+
     ylab("log of mean expression in activated cell")+
     xlab("Number of lists genes are found in")+
     theme_light()+
     scale_fill_discrete("RIGs lists")
ggsave(Fig2C, file="Figure2C.pdf", width=7, height = 5)

limits <- quantile(genidf[gene_biotype=="protein_coding"][i.Mean>0, i.Mean], c(0.1, 0.9))
genidf[, ExpressionGroups:=factor(ifelse(i.Mean<=0,"Not expressed", 
                                  ifelse(i.Mean<=limits[1], "Low10%",
                                         ifelse(i.Mean<=limits[2], "Mid","High10%")
                                         )
                                  ),
                                   levels=c("Not expressed", "Low10%", "Mid", "High10%")     
                                   )]
Figure2D <- ggplot(genidf[gene_biotype=="protein_coding"], aes(x=RIGSLISTE,fill=SuperenhancersA))+
     geom_bar(position="fill")+
     facet_grid(~ExpressionGroups)+
     theme_light()+
     theme(legend.position = "bottom")+
     ylab("")+
     xlab("")
ggsave(Figure2D, file="Figure2D.pdf", width=10, height = 5)

Fig3D <- ggplot(gg[gg$gene_biotype=="protein_coding",], aes(numberOfLists,Mean))+
     geom_violin()+
     scale_y_continuous()+
     ylab("Rlog of mean expression value")+
     theme_light()+
     xlab("Activated JQ1 cells")
ggsave(Fig3D, file="Figure3D.pdf", width=7, height = 5)


genidf$negLog <- -1*log(genidf$Activated1padjuponJQ1,base = 10)
SupplementaryFig2B <- ggplot(genidf[gene_biotype=="protein_coding"], aes(Activated1LFCuponJQ1, negLog, color=SuperenhancersA))+
     geom_point(alpha=0.5)+theme_light()+
     coord_cartesian(ylim=c(0,80))+
     ylab("-log10FC(p adjusted value)")+
     xlab("log2FC")+
     theme(legend.position = "bottom")
ggsave(SupplementaryFig2B, file="SupplementaryFigure2B.pdf", width=7, height = 7)


library(RColorBrewer)
genidf$ChangeInActivatedUponJQ1
SupplementaryFig2C <- ggplot(genidf[gene_biotype=="protein_coding"], aes(RIGSLISTE, fill=ChangeInActivatedUponJQ1))+
     geom_bar(position="fill")+
     facet_wrap(~SuperenhancersA)+
     ylab("Percentage of genes")+
     xlab("")+
     scale_fill_manual(
          name = "\nBehavior upon JQ1 treatment",
          values = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(3)[c(1,2,3)],
          breaks = c("DownregulatedUponJQ1","UnchangedUponJQ1","UpregulatedUponJQ1"),
          labels = c("Downregulated", "No significant change", "Upregulated")
     )
ggsave(SupplementaryFig2C, file="SupplementaryFigure2C.pdf", width=7, height = 5)
