#!/usr/bin/env Rscript

# shared_ancestry_simulator.R
# Generate models using ms and seq-gen based on model specified by YAML file

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# October-December 2013

library(tools)
library(parallel)
library(plyr)
library(yaml)
library(optparse)
library(ggplot2,quietly=TRUE)
library(gridExtra,quietly=TRUE)
suppressMessages(library(reshape,quietly=TRUE))

# phyclust contains ms and seq-gen,
# but also depends on ape,
# making the DNAbin class available
library(phyclust, quietly=TRUE)

# Parse options
options<-list(make_option(c("-T","--threads"),default=1,help="Number of threads [default %default]",metavar="numthreads"),
              make_option(c("-w","--windows"),default=100,help="Number of windows to simulate",metavar="numwindows"),
              make_option(c("-c","--combined"),type="character",help="YAML file list with proportions for combined models, eg -c model1.yml:0.1,model2.yml:0.9")
         )

opt<-parse_args(OptionParser(option_list=options))

if (is.null(opt$combined)) stop("Please specify model files in YAML format and proportions with -c, eg -c model1.yml:0.1,model2.yml:0.9")
if(!is.numeric(opt$windows) | opt$windows<1) stop("Please specify a positive number of windows with -w")
if(!is.numeric(opt$threads) | opt$threads<1) stop("Please specify a positive number of threads with -T")

#Requires DNAbin input and vector of population names for each sequence in the sample
nucdiv <- function(aln, popIndex) {
  #check if pops were defined
  if (is.null(popIndex) == TRUE) {
    stop("")
    pops <- list( "1" = labels(aln))
    } else {
    pops <- list()
    for (popName in unique(popIndex)) {
      pops[[popName]] <- labels(aln)[which(popIndex == popName)]
    }
  }
  #create an output matrix. It will have all populations as rows and columns, so upper and lower triangles will be mirrored.
  #this just makes indexing easier. Each block will be dxy for that pair pop populations, except diagnols which wll be pi.
  output <- matrix(nrow = length(pops), ncol = length (pops))
  rownames(output) <- names(pops)
  colnames(output) <- names(pops)
  #Get pairwise distance matrix 
  pdm <- dist.dna(aln, model = "raw", pairwise.deletion=TRUE, as.matrix = TRUE)
  #Now, for each pop, get pi and for each pair get dxy
  #Pi is the mean of the lower triangle, exclusing the diagnol, when the same pop is used to index the rows and columns:
  for (X in 1:length(pops)){
    for (Y in X:length(pops)){
      sub_pdm <- pdm[pops[[X]],pops[[Y]]]
      #remove upper triangle and, if they're the same pop, remove diagnol, then take the mean and add to output matrix
      if (X == Y){
        sub_pdm[upper.tri(sub_pdm, diag = TRUE)] <- NA
        output[X,Y] <- mean(sub_pdm, na.rm = T)
      } else {
        sub_pdm[upper.tri(sub_pdm, diag = FALSE)] <- NA
        output[X,Y] <- mean(sub_pdm, na.rm = T)
        output[Y,X] <- mean(sub_pdm, na.rm = T) # add to both lower and upper - for ease of indexing
      }
    }
  }
  output
}


nucdiv4pop<-function(aln,popIndex) {
    if (length(unique(popIndex)) != 4) stop("Please supply four populations only for nucleotide divergence calculations")
    ndmat<-nucdiv(aln,popIndex)
    data.frame(
        P1_Pi=ndmat[1,1],P2_Pi=ndmat[2,2],P3_Pi=ndmat[3,3],P4_Pi=ndmat[4,4],
        P1P2_dxy=ndmat[1,2],P1P3_dxy=ndmat[1,3],P1P4_dxy=ndmat[1,4],P2P3_dxy=ndmat[2,3],P2P4_dxy=ndmat[2,4],P3P4_dxy=ndmat[3,4]
    )
}


get.pop.slice<-function(popsites,bi.sites) {
    counts<-apply(popsites,2,base.freq,freq=TRUE)
    alleles<-apply(counts/counts,2,sum,na.rm=TRUE)
    s<-length(which(alleles>1))
    list(counts=counts,counts.bi=counts[,bi.sites],alleles=alleles,s=s)
}

get.abba.baba.stats<-function(counts,bi,P1,P2,P3,P4) {
    ABBA = 0
    BABA = 0
    maxABBA = 0
    maxBABA = 0
    maxABBA_B = 0
    maxBABA_B = 0

    for (i in 1:length(bi)){
      i.bi.counts<-counts[,bi[i]]
      alleles <- names(which(i.bi.counts > 0))
      anc <- names(which(P4$counts.bi[,i] != 0 ))
      if (length(anc) != 1) { anc <- names(which(i.bi.counts == max(i.bi.counts)))[1] }
      derived = alleles[which(alleles != anc)]
      P1df <- P1$counts.bi[derived,i] / sum(P1$counts.bi[,i])
      P2df <- P2$counts.bi[derived,i] / sum(P2$counts.bi[,i])
      P3df <- P3$counts.bi[derived,i] / sum(P3$counts.bi[,i])
      P4df <- P4$counts.bi[derived,i] / sum(P4$counts.bi[,i])
      ABBA <- sum(ABBA, (1 - P1df) * P2df * P3df * (1 - P4df), na.rm = TRUE)
      BABA <- sum(BABA, P1df * (1 - P2df) * P3df * (1 - P4df), na.rm = TRUE)
      maxABBA <- sum(maxABBA, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
      maxBABA <- sum(maxBABA, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
      if (is.na(P3df) == FALSE & is.na(P2df) == FALSE & P3df >= P2df){
        maxABBA_B <- sum(maxABBA_B, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
        maxBABA_B <- sum(maxBABA_B, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
        } else{
        maxABBA_B <- sum(maxABBA_B, (1 - P1df) * P2df * P2df * (1 - P4df), na.rm = TRUE)
        maxBABA_B <- sum(maxBABA_B, P1df * (1 - P2df) * P2df * (1 - P4df), na.rm = TRUE)
      }
    }
    output<-data.frame(ABBA=ABBA,BABA=BABA,
               maxABBA=maxABBA, maxBABA=maxBABA,
               maxABBA_B=maxABBA_B, maxBABA_B=maxBABA_B,
               D=(ABBA - BABA) / (ABBA + BABA),
               f=(ABBA - BABA) / (maxABBA - maxBABA),
               mf=(ABBA - BABA) / (maxABBA_B - maxBABA_B)
    )
    output$mfD0<-apply(output,1,function(x) if ((is.na(x["D"])) || (as.numeric(x["D"]) < 0)) NA else as.numeric(x["mf"]))
    output
}

# Requires DNAbin input, vector of population names for each sequence in the sample,
# and names for P1, P2, P3 and P4

ABBABABA <- function(aln, popIndex, P1, P2, P3, P4) {

  #get allele counts to find variable sites
  counts.all<-apply(aln,2,base.freq,freq=TRUE)
  #get numbers of alleles at each site (counts/counts converts allele counts to 0 or 1)
  num.alleles.all<-colSums(counts.all/counts.all,na.rm=TRUE)

  #biallelic sites
  bi.all <- which(num.alleles.all == 2)

  #varying sites
  var.sites.all<-num.alleles.all>1
  aln.var<-aln[,var.sites.all]

  alleles.var<-num.alleles.all[var.sites.all]
  s.all<-length(alleles.var)

  bi.var<-which(alleles.var==2)

  P1.slice<-get.pop.slice(aln.var[popIndex==P1,],bi.var)
  P2.slice<-get.pop.slice(aln.var[popIndex==P2,],bi.var)
  P3.slice<-get.pop.slice(aln.var[popIndex==P3,],bi.var)
  P4.slice<-get.pop.slice(aln.var[popIndex==P4,],bi.var)

  abba.baba.stats<-get.abba.baba.stats(counts.all,bi.all,P1.slice,P2.slice,P3.slice,P4.slice)

  data.frame(S=s.all, P1_S=P1.slice$s, P2_S=P2.slice$s, P3_S=P3.slice$s, P4_S=P4.slice$s, abba.baba.stats)
}

check.pops<-function(pops,totalpops,popnum) {
    if (length(pops) != popnum) stop(paste("Expected",popnum,"pops, got ",length(pops)))
    for (i in pops) {
        if (class(i) != "integer") stop("All population values must be integers")
        if (i > totalpops) stop("Population value larger than number of populations")
    }
}

parse.model<-function(filename) {
    model<-yaml.load_file(filename)

    if (!("pops" %in% names(model$ms))) stop("Please specify ms populations using the name 'pops'")
    if (!("opts" %in% names(model$ms))) stop("Please specify ms options using the name 'opts'")

    # Parse populations
    model$ms$popnum<-length(model$ms$pops)
    model$ms$nsam<-0
    model$ms$popnames<-vector()
    model$ms$popsizes<-rep(NA,model$ms$popnum)
    for (pop in 1:model$ms$popnum) {
        model$ms$nsam<-model$ms$nsam + model$ms$pops[[pop]]$seqs
        model$ms$popnames<-c(model$ms$popnames,rep(model$ms$pops[[pop]]$name,model$ms$pops[[pop]]$seqs))
        model$ms$popsizes[pop]<-model$ms$pops[[pop]]$seqs
    }
    
    # Validate migrations
    if ("migration" %in% names(model$ms$opts)) {
        if (!("n0" %in% names(model$ms$opts))) stop("Found migration options but no N0; please specify n0 ms option")
        for (i in model$ms$opts$migration) {
            if (!("mij" %in% names(i))) stop("Please specify mij for all migration options")
            if (class(i$mij) != "numeric") stop("migration mij parameters should be numeric")

            if (!("pops" %in% names(i))) stop ("Please specify pops for all migration options")
            check.pops(i$pops,length(model$ms$pops),2)
        }
    }

    # Validate splits
    if ("split" %in% names(model$ms$opts)) {
        if (!("n0" %in% names(model$ms$opts))) stop("Found split options but no N0; please specify n0 ms option")
        if (!("genperyear" %in% names(model$ms$opts))) stop("Found split options but no generations per year; please specify genperyear ms option")
        for (i in model$ms$opts$split) {
            if (!("time" %in% names(i)) & !("units" %in% names(i))) stop("Please specify time or units for all split options")
            if (!("pops" %in% names(i))) stop("Please specify pops for all split options")
            check.pops(i$pops,length(model$ms$pops),2)
        }
    }

    # Validate theta
    if (!("s" %in% names(model$seqgen))) {
        if (!("n0" %in% names(model$ms$opts))) stop("No seqgen s and no n0 provided for calculation; please specify seqgen s or ms n0 and ms mu")
        if (!("mu" %in% names(model$ms$opts))) stop("No seqgen s and no mu provided for calculation; please specify seqgen s or ms n0 and ms mu")
    }

    # Validate N0
    if ("n0" %in% names(model$ms$opts)) {
        n0test<-eval(parse(text=model$ms$opts$n0))
        if (class(n0test) != "numeric") stop("n0 must be a number or an expression producing a number")
        if (length(n0test) != 1) stop("n0 should only be a single number")
    }

    model$ms$thetaopt<-""
    if ("theta" %in% names(model$ms$opts)) {
        model$ms$thetaopt<-paste("-t",model$ms$opts$theta)
    }

    model$ms$sopt<-""
    if ("s" %in% names(model$ms$opts)) {
        model$ms$sopt<-paste("-s",model$ms$opts$s)
    }
    
    model$ms$optstring<-paste(c("-I",model$ms$popnum,model$ms$popsizes,"-T"),collapse=" ")
    if (model$ms$thetaopt != "") model$ms$optstring<-paste(model$ms$optstring,model$ms$thetaopt)
    if (model$ms$sopt != "") model$ms$optstring<-paste(model$ms$optstring,model$ms$sopt)
    model
}

get.model.options<-function(model) {
    
    model$ms$opts$n0val<-eval(parse(text=model$ms$opts$n0))

    if ("migration" %in% names(model$ms$opts)) {
        for (i in model$ms$opts$migration) {
            mig.rate<-i$mij * 4 * model$ms$opts$n0val
            migopt<-paste(c("-m",i$pops,mig.rate),collapse=" ")
            model$ms$optstring<-paste(model$ms$optstring,migopt)
        }
    }

    if ("split" %in% names(model$ms$opts)) {
        for (i in model$ms$opts$split) {
            if ("units" %in% names(i)) {
                split.time<-i$units
            } else {
                split.time<-(i$time * 1000000 * model$ms$opts$genperyear) / (4 * model$ms$opts$n0val)
            }
            splitopt<-paste(c("-ej",split.time,i$pops),collapse=" ")
            model$ms$optstring<-paste(model$ms$optstring,splitopt)
        }
    }

    if (!("s" %in% names(model$seqgen))) {
        model$seqgen$s<-4 * model$ms$opts$n0val * model$ms$opts$mu
    }

    model
}

simulate.window<-function(win,model) {
    model<-get.model.options(model)

    mstree<-ms(model$ms$nsam,1,model$ms$optstring)[3]
    simseqs<-seqgen(opts=sprintf("-m%s -l %d -s %f -q",model$seqgen$m,model$seqgen$l,model$seqgen$s),newick.tree=mstree)[-1]

    # Convert Seq-Gen output to DNAbin object
    bases<-matrix(nrow=model$ms$nsam,ncol=model$seqgen$l,byrow=TRUE)
    rownames(bases)<-1:model$ms$nsam
    for (s in 1:model$ms$nsam) {
        seqstats<-unlist(strsplit(simseqs[s],"\\s+"))
        seqstats[1]<-gsub('s','',seqstats[1])
        bases[as.numeric(seqstats[1]),]<-unlist(strsplit(tolower(seqstats[2]),''))
    }
    bases<-as.DNAbin(bases)
    
    windf<-data.frame(blockNumber=win,blockName=win,sites=model$seqgen$l)
    pardf<-data.frame(N0=model$ms$opts$n0val,s=model$seqgen$s)
    nddf<-nucdiv4pop(bases,model$ms$popnames)
    abdf<-ABBABABA(bases,model$ms$popnames,1,2,3,4)
    cbind(windf,pardf,nddf,abdf)
}

simulate.model<-function(filename,proportion=1,cumprop=0) {
    propwin<-ceiling(opt$windows*proportion)
    startwin<-ceiling(opt$windows*cumprop)+1
    endwin<-startwin+propwin-1
    if (file.access(filename)==-1) stop(sprintf("Model file %s does not exist",filename))
    model<-parse.model(filename)
    rbind.fill(mclapply(startwin:endwin,function(win) simulate.window(win,model),mc.cores=opt$threads))
}

propstats<-NULL
if (!(is.null(opt$combined))) {
    prop<-unlist(strsplit(opt$combined,","))
    modelfiles<-sapply(prop,function(p) unlist(strsplit(p,":"))[1])
    modelprops<-as.numeric(sapply(prop,function(p) unlist(strsplit(p,":"))[2]))
    opt$propdf<-data.frame(Filename=modelfiles,Proportion=modelprops,CumProp=cumsum(modelprops)-modelprops,stringsAsFactors=FALSE)
    if (sum(opt$propdf$Proportion) != 1) stop("Combined model proportions must add up to 1")

    propstats<-merge_all(mapply(simulate.model, opt$propdf$Filename,opt$propdf$Proportion,opt$propdf$CumProp,SIMPLIFY=FALSE))

    propmodel<-unlist(mapply(function(f,p) rep(f,p*opt$windows),opt$propdf$Filename,opt$propdf$Proportion))
    propstats<-cbind(propstats,Model=as.character(propmodel))

    propfilestub<-paste0(c(file_path_sans_ext(basename(opt$propdf$Filename)),opt$windows),collapse=".")
    write.csv(propstats,paste0(propfilestub,".csv"))

}
warninglist<-warnings()
if (!(is.null(warninglist))) warninglist