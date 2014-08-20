#!/usr/bin/env Rscript

# generate_summary_statistics.R
# Generate dxy and partitioning statistics for models generated using run_model_combinations.py

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013, May 2014, August 2014

library(parallel)
suppressMessages(library(reshape))
library(plyr)
library(optparse)
library(stringr)
library(tools)

options<-list(make_option(c("-t","--threads"),default=1,help="Number of threads [default %default]",metavar="numthreads"),
              make_option(c("-m","--modelfiles"),default="",type="character",help="Folder of CSV files with model output"),
              make_option(c("-r","--realdata"),default="",type="character",help="CSV file containing windows from real genome"),
              make_option(c("-l","--listofmodels"),default="",type="character",help="CSV file of model parameter combinations")
         )

opt<-parse_args(OptionParser(option_list=options))

if (opt$listofmodels == "" & opt$realdata == "") stop("Please specify a list of model parameter combinations in CSV format with -l")
if (opt$modelfiles == "" & opt$realdata == "") stop("Please specify a folder of model files with -m or a CSV file containing real data with -r")
if (opt$modelfiles != "" & !file.exists(opt$modelfiles)) stop(paste("Directory", opt$modelfiles, "does not exist! Please specify a valid directory with -m"))
if (opt$realdata != "" & !file.exists(opt$realdata)) stop(paste("Real data file", opt$realdata, "does not exist! Please specify a valid file with -r"))

partition.x<-function(model,stat,threshold=0.1,win) {
    outnum<-threshold*win
    Partition<-sapply(rank(-model[,stat],ties.method="random"), function(x) {if (x<=outnum) "_Outlier" else "_Background"})
    Partition<-paste0(stat,Partition)
    model<-cbind(model,Partition)
    names(model)[names(model)=="Partition"]<-paste0(stat,"_Partition")
    model
}

run.tests<-function(values,groups,test.name,test.func) {
    # Assumes order of groups is Background, Outlier
    formula<-values~groups
    two.tailed<-test.func(formula,alternative="two.sided",conf.int=TRUE)
    if (length(two.tailed$estimate)==1) mean.diff<-two.tailed$estimate else mean.diff<-two.tailed$estimate[1] - two.tailed$estimate[2]
    out<-data.frame(
        mean.diff=mean.diff,
        p.two.tailed=two.tailed$p.value,
        p.background.higher=test.func(formula,alternative="greater")$p.value,
        p.background.lower=test.func(formula,alternative="less")$p.value
    )
    names(out)<-paste(test.name,names(out),sep=".")
    out
}

summarise.tests<-function(model.sum,model.stats,model.partition,test.name,test.func) {
    model.tests<-t(apply(model.stats,2,function(x) unlist(run.tests(x,model.partition,test.name,test.func))))
    model.tests<-cbind(Variable=rownames(model.tests),model.tests)
    merge(model.sum,model.tests,by="Variable",sort=FALSE)
}

get.model.summary<-function(models,column,pars=NULL) {
    model.stats<-models[c("P1P2_dxy","P1P3_dxy","P2P3_dxy")]
    model.partition<-models[,column]
    
    model.means<-melt(aggregate(model.stats,list(model.partition),mean,na.rm=TRUE),id="Group.1")
    model.sds<-melt(aggregate(model.stats,list(model.partition),sd,na.rm=TRUE),id="Group.1")
    model.ses<-melt(aggregate(model.stats,list(model.partition), function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}),id="Group.1")
    model.sum<-cbind(pars,model.means,model.sds$value,model.ses$value)
    names(model.sum)<-c(names(pars),"Model","Variable","Mean","SD","SE")

    model.sum<-summarise.tests(model.sum,model.stats,model.partition,"t",t.test)
    model.sum<-summarise.tests(model.sum,model.stats,model.partition,"w",wilcox.test)

    model.sum
}

calc.fD0<-function(df, fstat) {
    apply(df,1,function(x) if ((is.na(x["D"])) || (as.numeric(x["D"])<0)) NA else as.numeric(x[fstat]))
}

get.fs<-function(models, altprop, windows) {

    models$fGD0<-calc.fD0(models,"fG")
    models$fhomD0<-calc.fD0(models,"fhom")
    models$fdD0<-calc.fD0(models,"fd")

    models<-partition.x(models,"D",altprop,windows)
    models<-partition.x(models,"fGD0",altprop,windows)
    models<-partition.x(models,"fhomD0",altprop,windows)
    models<-partition.x(models,"fdD0",altprop,windows)
    
    models
}

mean.pc<-function(stat) {
    mean(stat, na.rm=TRUE) * 100
}

sd.pc<-function(stat) {
    sd(stat, na.rm=TRUE) * 100
}

summarise.model<-function(model, modeldir, seqs) {
    filestub = paste0("Alternate_t123-", model$Alternate_t123, "_t23-", model$Alternate_t23, ".Background_t123-", model$Background_t123, "_t21-", model$Background_t21, ".*.", seqs, ".csv")
    filename = dir(modeldir, filestub, full.names=TRUE)
    pars<-data.frame(str_match(filename,"Alternate_t123-(.+)_t23-(.+)\\.Background_t123-(.+)_t21-(.+)\\.(.+)\\..+\\.csv"),stringsAsFactors=FALSE)
    names(pars)<-c("File","Alternate_t123","Alternate_t23","Background_t123","Background_t21","Windows")
    pars$Windows<-as.numeric(pars$Windows)

    models<-read.csv(filename)
    models$Model<-factor(str_match(models$Model,"(.+)_t123")[,2],levels=c("Background","Alternate"))
    models$Model<-revalue(models$Model,c(Background="Real_Background",Alternate="Real_Outlier"))

    models<-get.fs(models, 0.1, pars$Windows)
    
    model.summary<-get.model.summary(models,"Model",pars)
    D.summary<-get.model.summary(models,"D_Partition",pars)
    fGD0.summary<-get.model.summary(models,"fGD0_Partition",pars)
    fhomD0.summary<-get.model.summary(models,"fhomD0_Partition",pars)
    fdD0.summary<-get.model.summary(models,"fdD0_Partition",pars)

    model.sum<-rbind(model.summary, D.summary, fGD0.summary, fhomD0.summary, fdD0.summary)

    is.D.positive<-models$D>=0
    real.outliers<-models$Model=="Real_Outlier"
    
    model.comp<-data.frame(pars,
        Stat=c("D","fGD0","fhomD0","fdD0"),
        Stat.Mean=c(mean.pc(models$D), mean.pc(models$fGD0), mean.pc(models$fhomD0), mean.pc(models$fdD0)),
        Stat.SD=c(sd.pc(models$D), sd.pc(models$fGD0), sd.pc(models$fhomD0), sd.pc(models$fdD0)),
        Dpos=sum(is.D.positive,na.rm=TRUE),
        DposAlternate=sum(is.D.positive & real.outliers,na.rm=TRUE),
        DposOutlier=c(sum(is.D.positive & models$D_Partition=="D_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$fGD0_Partition=="fGD0_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$fhomD0_Partition=="fhomD0_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$fdD0_Partition=="fdD0_Outlier", na.rm=TRUE)),
        OutlierAlternate=c(sum(real.outliers & models$D_Partition=="D_Outlier",na.rm=TRUE),
                           sum(real.outliers & models$fGD0_Partition=="fGD0_Outlier",na.rm=TRUE),
                           sum(real.outliers & models$fhomD0_Partition=="fhomD0_Outlier",na.rm=TRUE),
                            sum(real.outliers & models$fdD0_Partition=="fdD0_Outlier",na.rm=TRUE))
    )
    model.comp$OutlierAlternatePC<-model.comp$OutlierAlternate/model.comp$DposOutlier*100

    list(model.sum,model.comp)
}

summarise.null<-function(model, modeldir, seqs) {
    filepattern = paste0("^Background_t123-", model$Background_t123, "_t21-", model$Background_t21,".*.", seqs, ".csv")
    filename = dir(modeldir, filepattern, full.names=TRUE)
    pars<-data.frame(str_match(filename,"Background_t123-(.+)_t21-(.+)\\.(.+)\\..+\\.csv"),stringsAsFactors=FALSE)
    names(pars)<-c("File","Background_t123","Background_t21","Windows")
    pars$Windows<-as.numeric(pars$Windows)

    models<-read.csv(filename)

    models<-get.fs(models, 0.1, pars$Windows)
    
    D.summary<-get.model.summary(models,"D_Partition",pars)
    fGD0.summary<-get.model.summary(models,"fGD0_Partition",pars)
    fhomD0.summary<-get.model.summary(models,"fhomD0_Partition",pars)
    fdD0.summary<-get.model.summary(models,"fdD0_Partition",pars)

    model.sum<-rbind(D.summary, fGD0.summary, fhomD0.summary, fdD0.summary)

    is.D.positive<-models$D>=0
    model.comp<-data.frame(pars,
        Stat=c("D","fGD0","fhomD0","fdD0"),
        Stat.Mean=c(mean.pc(models$D), mean.pc(models$fGD0), mean.pc(models$fhomD0), mean.pc(models$fdD0)),
        Stat.SD=c(sd.pc(models$D), sd.pc(models$fGD0), sd.pc(models$fhomD0), sd.pc(models$fdD0)),
        Dpos=sum(is.D.positive,na.rm=TRUE),
        DposOutlier=c(sum(is.D.positive & models$D_Partition=="D_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$fGD0_Partition=="fGD0_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$fhomD0_Partition=="fhomD0_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$fdD0_Partition=="fdD0_Outlier", na.rm=TRUE))
    )

    list(model.sum,model.comp)
}

process.alt<-function(modeldir, model.df, seqs) {
    model.list<-split(model.df, 1:nrow(model.df)) # Process df by rows
    csv.summary<-mclapply(model.list, summarise.model, modeldir=modeldir, seqs=seqs, mc.cores=opt$threads)

    models<-rbind.fill(lapply(csv.summary,function(x) x[[1]]))
    model.comp<-rbind.fill(lapply(csv.summary,function(x) x[[2]]))

    models<-merge(models,unique(model.df[,c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","ModelType")]),by=c("Background_t123","Background_t21","Alternate_t123","Alternate_t23"))
    write.table(models,file=paste0(opt$modelfiles,".alternate_models.dxy.summary.", seqs, ".tsv"),sep="\t",row.names=FALSE)

    model.comp<-merge(model.comp,unique(model.df[,c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","ModelType")]),by=c("Background_t123","Background_t21","Alternate_t123","Alternate_t23"))
    write.table(model.comp,file=paste0(opt$modelfiles,".alternate_models.partition.summary.", seqs, ".tsv"),sep="\t",row.names=FALSE)    
}

process.null<-function(modeldir, model.df, seqs) {
    model.list<-split(model.df, 1:nrow(model.df)) # Process df by rows
    csv.summary<-mclapply(model.list, summarise.null, modeldir=modeldir, seqs=seqs, mc.cores=opt$threads)
    
    models<-rbind.fill(lapply(csv.summary,function(x) x[[1]]))
    model.comp<-rbind.fill(lapply(csv.summary,function(x) x[[2]]))

    models<-merge(models,unique(model.df[,c("Background_t123","Background_t21","ModelType")]),by=c("Background_t123","Background_t21"))
    write.table(models,file=paste0(opt$modelfiles,".null_models.dxy.summary.", seqs, ".tsv"),sep="\t",row.names=FALSE)

    model.comp<-merge(model.comp,unique(model.df[,c("Background_t123","Background_t21","ModelType")]),by=c("Background_t123","Background_t21"))
    write.table(model.comp,file=paste0(opt$modelfiles,".null_models.partition.summary.", seqs, ".tsv"),sep="\t",row.names=FALSE)
    
}

options(width=200)

if (opt$modelfiles != "") {
    
    model.list<-read.csv(opt$listofmodels,stringsAsFactors=FALSE,colClasses="character")
    
    for (seqs in c("ms","sg")) {
        process.alt(opt$modelfiles, model.list[model.list$ModelType != "Null model",], seqs)
        process.null(opt$modelfiles, model.list[model.list$ModelType == "Null model",], seqs)
    }
}

if (opt$realdata != "") {
    real<-read.csv(opt$realdata)
    real$fdD0<-calc.fD0(real, "fd")
    real<-partition.x(real,"D",0.1,nrow(real))
    real<-partition.x(real,"fGD0",0.1,nrow(real))
    real<-partition.x(real,"fhomD0",0.1,nrow(real))
    real<-partition.x(real,"fdD0",0.1,nrow(real))
    
    pars<-data.frame(File=opt$realdata,Windows=nrow(real))
    real.stats<-rbind(get.model.summary(real,"D_Partition",pars),get.model.summary(real,"fGD0_Partition",pars),get.model.summary(real,"fhomD0_Partition",pars),get.model.summary(real,"fdD0_Partition",pars))
    write.table(real.stats,file=paste0(file_path_sans_ext(basename(opt$realdata)),".dxy.summary.tsv"),sep="\t",row.names=FALSE)
}

warninglist<-warnings()
if(!(is.null(warninglist))) warninglist