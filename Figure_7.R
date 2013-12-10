#!/usr/bin/env Rscript

library(optparse)
suppressMessages(library(ggplot2))
library(plyr)

options(width=1000)

options<-list(
    make_option(c("-i","--input"),default=NULL,type="character",help="TSV file with model Dxy summary statistics")
)

opt<-parse_args(OptionParser(option_list=options))
if (is.null(opt$input)) stop("Please supply an input file with -i")

match.pred.stat<-function(models,stat) {
    stats<-c(paste0(stat,"_Outlier"),paste0(stat,"_Background"))
    mdxy.stat<-subset(models,Model %in% stats)
    mdxy.p2p3<-mdxy.stat[mdxy.stat$Variable=="P2P3_dxy",]
    data.frame(Stat=stat,
               P2P3_dxy.Outlier   =mdxy.p2p3$Mean[2],
               P2P3_dxy.Background=mdxy.p2p3$Mean[1],
               P2P3_dxy.effect    =mdxy.p2p3$Mean[2]/mdxy.p2p3$Mean[1]*100
              )
}

get.dxy.props<-function(models) {
    mdxy<-subset(models,grepl("dxy",Variable))
    ab<-rbind(
        match.pred.stat(mdxy,"Real"),
        match.pred.stat(mdxy,"D"),
        match.pred.stat(mdxy,"mfD0")
    )
    ab$Stat<-revalue(factor(ab$Stat,levels=c("Real","D","mfD0")),c(Real="Simulation",mfD0="f"))
    ab
}

models<-read.delim(opt$input,stringsAsFactors=FALSE)

stats<-ddply(models,.(Background_t123,Background_t21,Alternate_t123,Alternate_t23,ModelType),get.dxy.props)

stats$ModelType<-revalue(stats$ModelType,c("Ancestral structure"="Ancestral Structure","gene flow P2-P3"="Gene Flow P2-P3","gene flow P3-P2"="Gene Flow P3-P2"))

stats$TreeDistProp<-apply(stats,1,function(x) {
                params<-as.numeric(c(x["Background_t123"],x["Background_t21"],x["Alternate_t123"],x["Alternate_t23"]))
                (max(params) - min(params))/max(params)
            })


d.f.dxy<-by(stats,stats$ModelType,function(x) cbind(x$ModelType,x["Background_t123"],x["Background_t21"],x["Alternate_t123"],x["Alternate_t23"],prop=x[x$Stat=="D",]$P2P3_dxy.effect-x[x$Stat=="f",]$P2P3_dxy.effect))

pdf("Figure_7.pdf",width=3.25,height=4.5)
ggplot(stats,aes(ModelType,P2P3_dxy.effect,colour=ModelType))+
    geom_point(size=1.3)+
    facet_grid(.~Stat)+
    scale_y_continuous(limits=c(0,110),breaks=c(0,20,40,60,80,100))+
    guides(colour=FALSE)+
    labs(x="Model Type", y=expression(paste("Outlier to Non-Outlier P2-P3 ",italic(d[XY])," (%)")), colour="Model Type")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text=element_text(size=10),
          axis.title.y=element_text(size=10))

dev.off()

warninglist<-warnings()
if (!(is.null(warninglist))) warninglist