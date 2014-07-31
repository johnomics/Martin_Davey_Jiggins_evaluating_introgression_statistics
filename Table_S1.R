#!/usr/bin/env Rscript

# Table_S1.R
# Generate table S1 of all model results from summary statistics

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013, May 2014

library(plyr)

options(width=1000, stringsAsFactors=FALSE)

match.dxy.stat<-function(models,stat) {
    stats<-c(paste0(stat,"_Outlier"),paste0(stat,"_Background"))
    mdxy.stat<-subset(models,Model %in% stats)
    mdxy.p2p3<-mdxy.stat[mdxy.stat$Variable=="P2P3_dxy",]
    data.frame(Stat=stat,
               P2P3_dxy.Outlier.Mean    = round(mdxy.p2p3$Mean[2], 7),
               P2P3_dxy.Outlier.SD      = round(mdxy.p2p3$SD[2]  , 7),
               P2P3_dxy.Background.Mean = round(mdxy.p2p3$Mean[1], 7),
               P2P3_dxy.Background.SD   = round(mdxy.p2p3$SD[1]  , 7),
               P2P3_dxy.Significant     = mdxy.p2p3$w.p.background.higher[1]<(0.01/mdxy.p2p3$ModelCount[1])
              )
}

get.dxy.stats<-function(dxy.summary) {
    ab<-rbind(
        match.dxy.stat(dxy.summary,"Real"),
        match.dxy.stat(dxy.summary,"D"),
        match.dxy.stat(dxy.summary,"fGD0"),
        match.dxy.stat(dxy.summary,"fhomD0"),
        match.dxy.stat(dxy.summary,"fdD0")
    )
    ab$Stat<-revalue(factor(ab$Stat,levels=c("Real","D","fGD0","fhomD0","fdD0")),c(Real="Simulation",fGD0="fG", fhomD0="fhom", fdD0="fd"))
    ab
}

get.partition.stat<-function(partition.summary) {
    data.frame(Stat=c("Simulation","D", "fG", "fhom", "fd"),
               Outlier.Alternate.PC = c(
                   NA,
                   partition.summary[partition.summary$Stat=="D",]$OutlierAlternatePC,
                   partition.summary[partition.summary$Stat=="fGD0",]$OutlierAlternatePC,
                   partition.summary[partition.summary$Stat=="fhomD0",]$OutlierAlternatePC,
                   partition.summary[partition.summary$Stat=="fdD0",]$OutlierAlternatePC
               )
    )
}

get.model.summary<-function(fileprefix, filepostfix, recombval) {
    alt.dxy.summary.file<-paste0(fileprefix, "alternate_models.dxy", filepostfix)
    alt.partition.summary.file<-paste0(fileprefix, "alternate_models.partition", filepostfix)
    null.dxy.summary.file<-paste0(fileprefix, "null_models.dxy", filepostfix)
    null.partition.summary.file<-paste0(fileprefix, "null_models.partition", filepostfix)

    alt.dxy.summary<-read.delim(alt.dxy.summary.file)
    null.dxy.summary<-read.delim(null.dxy.summary.file)
    null.dxy.summary<-cbind(null.dxy.summary[1:2],Alternate_t123=NA, Alternate_t23=NA,null.dxy.summary[3:length(null.dxy.summary)])
    dxy.summary<-rbind(alt.dxy.summary, null.dxy.summary)


    alt.partition.summary<-read.delim(alt.partition.summary.file)
    null.partition.summary<-read.delim(null.partition.summary.file)
    null.partition.summary<-cbind(null.partition.summary[1:2],
                                  Alternate_t123=NA, Alternate_t23=NA,
                                  null.partition.summary[3:8],
                                  DposAlternate=NA,
                                  null.partition.summary[9],
                                  OutlierAlternate=NA, OutlierAlternatePC=NA,
                                  null.partition.summary[10])
    partition.summary<-rbind(alt.partition.summary, null.partition.summary)

    model.counts<-table(unique(subset(dxy.summary,select=c("File","ModelType")))$ModelType)

    dxy.summary$ModelCount<-sapply(dxy.summary$ModelType, function(x) model.counts[x])
    dxy.stats<-ddply(dxy.summary,
                     .(Background_t123, Background_t21, Alternate_t123, Alternate_t23, ModelType, ModelCount),
                     get.dxy.stats
                    )
    partition.stats<-ddply(partition.summary,
                           .(Background_t123, Background_t21, Alternate_t123, Alternate_t23, ModelType),
                           get.partition.stat
                    )
    stats<-merge(dxy.stats,partition.stats)

    stats$ModelCount<-NULL
    names(stats)[names(stats)=="Background_t123"]<-"t23"
    names(stats)[names(stats)=="Background_t21"]<-"t12"
    stats<-cbind(stats[1:2],
                 tGF  = ifelse(grepl("gene flow", stats$ModelType), stats$Alternate_t23, NA),
                 tSTR = ifelse(stats$ModelType=="Ancestral structure", stats$Alternate_t123, NA),
                 stats[5:length(stats)]
                )

}

r50<-cbind("4Nr"=0.01, get.model.summary("model_files_win10000_s0.01_l5000_r50.", ".summary.sg.tsv", 50))
r5<-cbind("4Nr"=0.001,get.model.summary("model_files_win10000_s0.01_l5000_r5.", ".summary.sg.tsv", 5))

write.table(rbind(r50,r5), file="Table_S1.tsv", sep="\t", quote=FALSE, row.names=FALSE)