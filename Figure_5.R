#!/usr/bin/env Rscript

# Figure_5.R
# Generate the plot for Figure 5 as a PDF file

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013, May 2014

library(optparse)
suppressMessages(library(ggplot2))
library(plyr)

options(width=1000)

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
        match.pred.stat(mdxy,"fGD0"),
        match.pred.stat(mdxy,"fhomD0"),
        match.pred.stat(mdxy,"fdD0")
    )
    ab$Stat<-revalue(factor(ab$Stat,levels=c("Real","D","fGD0","fhomD0","fdD0")),c(Real="Simulation",fGD0="fG", fhomD0="fhom", fdD0="fd"))
    ab
}


make.stats<-function(altfile,nullfile, recombval) {
    alt.models<-read.delim(altfile,stringsAsFactors=FALSE)
    null.models<-read.delim(nullfile,stringsAsFactors=FALSE)

    null.models<-cbind(null.models[1:2],Alternate_t123=NA, Alternate_t23=NA,null.models[3:length(null.models)])

    models<-rbind(alt.models,null.models)

    stats<-ddply(models,.(Background_t123,Background_t21,Alternate_t123,Alternate_t23,ModelType),get.dxy.props)
    stats<-stats[!is.na(stats$P2P3_dxy.effect),]
    stats<-cbind(stats,Recombination=paste0("4Nr=",recombval))

    stats$ModelType<-revalue(factor(stats$ModelType),c("Null model"="Null Model","Ancestral structure"="Ancestral Structure","gene flow P2-P3"="Gene Flow P2-P3","gene flow P3-P2"="Gene Flow P3-P2"))
    stats$ModelType<-factor(stats$ModelType, levels=c("Gene Flow P2-P3", "Gene Flow P3-P2", "Ancestral Structure", "Null Model"))

    stats
}


r50.stats<-make.stats(altfile  = "model_files_win10000_s0.01_l5000_r50.alternate_models.dxy.summary.sg.tsv",
                      nullfile = "model_files_win10000_s0.01_l5000_r50.null_models.dxy.summary.sg.tsv",
                      recombval = 0.01)


r5.stats<-make.stats(altfile  = "model_files_win10000_s0.01_l5000_r5.alternate_models.dxy.summary.sg.tsv",
                      nullfile = "model_files_win10000_s0.01_l5000_r5.null_models.dxy.summary.sg.tsv",
                      recombval = 0.001)

stats<-rbind(r50.stats, r5.stats)

pdf("Figure_5.pdf",width=6.75,height=4.5)
ggplot(stats, aes(Stat, P2P3_dxy.effect, colour=ModelType)) +
    geom_point(size=1) +
    facet_grid(.~Recombination+ModelType, scales="free_x", space="free_x") +
    scale_x_discrete(labels=c("Simulation"="Simulation", "D"=expression(italic("D")),"fG"=expression(italic(f["G"])), "fhom"=expression(italic(f["hom"])), "fd"=expression(italic(f["d"])))) +
    scale_y_continuous(limits=c(0,110),breaks=c(0,20,40,60,80,100)) +
    scale_colour_manual(values=c("#00ba38", "#619cff", "#f8766d", "#ad5c08")) +
    guides(colour=FALSE) +
    labs(y=expression(paste("Outlier to Non-Outlier P2-P3 ",italic(d[XY])," (%)")), colour="Model Type") +
    theme_bw() +
    theme(axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text=element_text(size=6),
          axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
          axis.title.y=element_text(size=10),
          legend.position="none")

dev.off()

warninglist<-warnings()
if (!(is.null(warninglist))) warninglist