library(xlsx)
library(ggplot2)
TableM1 <- read.xlsx("Additional file 1.xlsx",1, header = T)
dd <- as.data.frame(t(TableM1[,2:9]))
nms <- c("Lucic", "Maldarelli", "Cohn", "Han","Ikeda", "Wagner", "Kok", "Brady")
dd <- dd[nms,]
dd$Source <- c("In vitro", "cART Patients", "cART Patients", "cART Patients", "cART Patients", "cART Patients", "cART Patients","In vitro")
ddd <- data.frame(pp=dd$V3, Source=dd$Source)

supplementary1a <- ggplot(ddd, aes(y=pp, x=Source, color=Source))+
    stat_boxplot(geom="errorbar", coef = 1.5)+
    geom_boxplot()+
    ylim(0,1)+
    theme_minimal()+
    ylab("% of unique IS in genes")+xlab("")
ggsave(supplementary1a, file="SupplementaryFigure1A.pdf", width = 7, height = 5)
