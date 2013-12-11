#!/usr/bin/env Rscript

# Figure_S3_S4.R
# Generate the plots for Figures S3 and S4 as PDF files

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013

suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
library(plyr)
library(stringr)

options(width=200)
load.with.window.size<-function(filename) {
    t<-read.delim(filename)
    t$Size<-as.numeric(str_match(filename,"size([[:digit:]]+)")[1,2])
    t[c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","File","Stat","ModelType","Size","OutlierAlternatePC")]
}

models<-rbind(load.with.window.size("model_files_win10000_size5000.partition.summary.tsv"),
              load.with.window.size("model_files_win5000_size10000.partition.summary.tsv"),
              load.with.window.size("model_files_win2500_size20000.partition.summary.tsv")
              )

kbwin<-function(x) {
    xkb<-as.numeric(x)/1000
    out<-paste(xkb,"kb windows")
    out
}

models$Stat<-revalue(models$Stat,c("mfD0"="f"))

models$Size<-factor(models$Size,levels=sort(unique(models$Size)))
models$Size<-mapvalues(models$Size,from=levels(models$Size),to=kbwin(levels(models$Size)))



models.f<-rename(subset(models,Stat=="f"),c("OutlierAlternatePC"="f_Accuracy"))
models.D<-rename(subset(models,Stat=="D"),c("OutlierAlternatePC"="D_Accuracy"))

models.figS3<-merge(models.f[c("File","ModelType","Size","f_Accuracy")],models.D[c("File","D_Accuracy")],by="File")

pdf("Figure_S3.pdf",10.5,3.5)

ggplot(models.figS3,aes(f_Accuracy,D_Accuracy,colour=ModelType))+
       facet_grid(.~Size)+
       geom_point()+
       geom_abline(intercept=0,slope=1,linetype="dashed")+
       xlim(c(0,100))+
       ylim(c(0,100))+
       xlab("% predicted by f")+
       ylab("% predicted by D")+
       coord_fixed()+
       theme_bw()
dev.off()

models.5kb<-rename(subset(models,Size=="5 kb windows"),c("OutlierAlternatePC"="Accuracy_5kb"))
models.10_20kb<-rename(subset(models,Size != "5 kb windows"),c("OutlierAlternatePC"="Accuracy_others"))

models.figS4<-merge(models.5kb[c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","Stat","ModelType","Accuracy_5kb")],
                    models.10_20kb[c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","Stat","Size","Accuracy_others")],
                    by=c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","Stat"))

pdf("Figure_S4.pdf",8,7)
ggplot(models.figS4,aes(Accuracy_5kb,Accuracy_others,colour=ModelType))+
       facet_grid(Size~Stat)+
       geom_point()+
       geom_abline(intercept=0,slope=1,linetype="dashed")+
       xlim(c(0,100))+
       ylim(c(0,100))+
       xlab("Accuracy (5 kb windows)")+
       ylab("Accuracy")+
       coord_fixed()+
       theme_bw()
dev.off()