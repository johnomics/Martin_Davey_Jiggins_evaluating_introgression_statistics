#!/usr/bin/env Rscript

# compare_f_estimators.R
# Simulations to compare D statistic and f estimators and make Figure 2 and Figure S1.

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# May 2014



library(phyclust)
library(parallel)
library(plyr)

ms.to.DNAbin <- function(ms_output, nsam=0, inc_inv=FALSE, len=NULL) {
    snp.strings <- ms_output[(length(ms_output)-nsam+1):length(ms_output)]
    binary.mat.var <- t(matrix(as.numeric(unlist(sapply(snp.strings,strsplit, split = ""))),ncol = length(snp.strings)))
    if (inc_inv==TRUE) {
        stopifnot(mode(len)=="numeric" & len >= length(binary.mat.var[1,]))
        nsam <- length(binary.mat.var[,1])
        nvar <- length(binary.mat.var[1,])
        binary.mat <- matrix(numeric(length=nsam*len), nrow=nsam)
        positions <- floor(as.numeric(unlist(strsplit(ms_output[length(ms_output)-nsam],split="    "))[2:(nvar+1)])*len)+1
        binary.mat[1:nsam,positions] <- binary.mat.var
    } else {
        binary.mat <- binary.mat.var
    }
    bases <- ifelse(binary.mat == "1", "t", "a")
    rownames(bases)<-1:nsam
    as.DNAbin(bases)
}



seqgen.to.DNAbin<-function(simseqs, len ) {

    bases<-matrix(nrow=length(simseqs),ncol=len,byrow=TRUE)
    rownames(bases)<-1:length(simseqs)
    for (s in 1:length(simseqs)) {
        seqstats<-unlist(strsplit(simseqs[s],"\\s+"))
        seqstats[1]<-gsub('s','',seqstats[1])
        bases[as.numeric(seqstats[1]),]<-unlist(strsplit(tolower(seqstats[2]),''))
    }
    as.DNAbin(bases)
}


get.pop.slice<-function(popsites,bi.sites) {
    counts<-apply(popsites,2,base.freq,freq=TRUE)
    alleles<-apply(counts/counts,2,sum,na.rm=TRUE)
    s<-length(which(alleles>1))
    list(counts=counts,counts.bi=counts[,bi.sites],alleles=alleles,s=s)
}



get.abba.baba.stats<-function(counts,bi,P1,P2,P3,P4,P31=NULL,P32=NULL) {
    
    do.real.f <- (is.null(P31)==F & is.null(P32)==F) #if both P3 subsets are given, calculate real f
    
    ABBA = 0
    BABA = 0
    maxABBA_G = 0
    maxBABA_G = 0
    maxABBA_hom = 0
    maxBABA_hom = 0
    maxABBA_D = 0
    maxBABA_D = 0

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
      if (do.real.f == T){
        P31df <- P31$counts.bi[derived,i] / sum(P31$counts.bi[,i])
        P32df <- P32$counts.bi[derived,i] / sum(P32$counts.bi[,i])
        }
      ABBA <- sum(ABBA, (1 - P1df) * P2df * P3df * (1 - P4df), na.rm = TRUE)
      BABA <- sum(BABA, P1df * (1 - P2df) * P3df * (1 - P4df), na.rm = TRUE)
      if (do.real.f==T){
        maxABBA_G <- sum(maxABBA_G, (1 - P1df) * P31df * P32df * (1 - P4df), na.rm = TRUE)
        maxBABA_G <- sum(maxBABA_G, P1df * (1 - P31df) * P32df * (1 - P4df), na.rm = TRUE)
          }
      maxABBA_hom <- sum(maxABBA_hom, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
      maxBABA_hom <- sum(maxBABA_hom, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
      if (is.na(P3df) == FALSE & is.na(P2df) == FALSE & P3df >= P2df){
        maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
        maxBABA_D <- sum(maxBABA_D, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
        } else{
        maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P2df * P2df * (1 - P4df), na.rm = TRUE)
        maxBABA_D <- sum(maxBABA_D, P1df * (1 - P2df) * P2df * (1 - P4df), na.rm = TRUE)
      }
    }
    output<-data.frame(ABBA=ABBA,BABA=BABA,
               maxABBA_G=maxABBA_G, maxBABA_G=maxBABA_G,
               maxABBA_D=maxABBA_D, maxBABA_D=maxBABA_D,
               D=(ABBA - BABA) / (ABBA + BABA),
               fG=(ABBA - BABA) / (maxABBA_G - maxBABA_G),
               fhom=(ABBA - BABA) / (maxABBA_hom - maxBABA_hom),
               fd=(ABBA - BABA) / (maxABBA_D - maxBABA_D)
               )
    output$fdD0<-apply(output,1,function(x) if ((is.na(x["D"])) || (as.numeric(x["D"]) < 0)) NA else as.numeric(x["fd"]))
    output
}



# Requires DNAbin input, vector of population names for each sequence in the sample,

# The arguments 'P1', 'P2', 'P3' and 'P4' (and optional pops 'P31' and 'P32') are
# each given as an integer vector of all sequences in each pop.
# eg. P1 = c(1,2,3,4,5,6,7,8), P2 = 9:16, P3 = 17:24, P4 = 25:32, P31 = 17:20, P32 = 21:24

ABBABABA <- function(aln, P1, P2, P3, P4, P31 = NULL, P32 = NULL) {
  

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

  P1.slice<-get.pop.slice(aln.var[P1,],bi.var)
  P2.slice<-get.pop.slice(aln.var[P2,],bi.var)
  P3.slice<-get.pop.slice(aln.var[P3,],bi.var)
  P4.slice<-get.pop.slice(aln.var[P4,],bi.var)
  if (is.null(P31)==F & is.null(P32)==F) {
    P31.slice<-get.pop.slice(aln.var[P31,],bi.var)
    P32.slice<-get.pop.slice(aln.var[P32,],bi.var)
    }

  abba.baba.stats<-get.abba.baba.stats(counts.all,bi.all,P1.slice,P2.slice,P3.slice,P4.slice,P31.slice,P32.slice)

  data.frame(S=s.all, P1_S=P1.slice$s, P2_S=P2.slice$s, P3_S=P3.slice$s, P4_S=P4.slice$s, abba.baba.stats)
}


sim_D_f <- function(n=1, nsam, P1, P2, P3, P4, P31, P32, ms_opts, len, scale, subModel, stats){
  ms.out<-ms(nsam=nsam, opts=ms_opts)
  ms.bases<-ms.to.DNAbin(ms.out, nsam)
  
  ms.trees<-ms.out[3:(grep("^segsites",ms.out)-1)]
  
  #simseqs<-system(sprintf("seq-gen -m%s -l %d -s %f -p %d -q", "HKY", len, scale, length(ms.trees)), input=ms.trees, intern=TRUE)[-1]
  
  seqgen.tmp.file<-paste("tmp", "seqgen", as.character(runif(1)), "out", sep=".")
  system(sprintf("seq-gen -m%s -l %d -s %f -p %d -q > %s", subModel, len, scale, length(ms.trees),  seqgen.tmp.file), input=ms.trees, ignore.stderr=TRUE)
  simseqs<-scan(seqgen.tmp.file, what="", sep="\n", quiet=TRUE)[-1]
  file.remove(seqgen.tmp.file)
  
  sg.bases<-seqgen.to.DNAbin(simseqs, len)
  
  abdf <- rbind(ABBABABA(ms.bases, P1,P2,P3,P4,P31,P32)[stats], ABBABABA(sg.bases, P1,P2,P3,P4,P31,P32)[stats])
  
  return(abdf[stats])
  }



#########################################################################################################################
### simulation. Multiple fs, multiple times of gene flow, using both ms and seqgen


subModel = "HKY"

P1 <- 1:8
P2 <- 9:16
P3 <- 17:24
P4 <- 25:32
P31 <- 17:20
P32 <- 21:24

nsam = 32

fs = seq(0,1.0,0.1)

nreps = 100

stats <- c("D","fG", "fhom","fd")

t = 100
len = 10000
rec_number = 100
tGF = c(0.1,0.5)
scale=0.01
theta = paste("-t", t)
pop_sizes = "-I 4 8 8 8 8"
join1 = "-ej 1 2 1"
join2 = "-ej 2 3 1"
join3 = "-ej 3 4 1"
rec = paste("-r", rec_number, len)


### gene flow 3-2
# three dimensional array for the output. Dimensions are fs (0.1,0.2 etc), stats (D, f etc), and reps (1:100).
all_data_32 <- array(dim = c(length(fs), nreps, length(stats), length(tGF), length(c("ms","seqgen"))), dimnames = list(fs, 1:nreps, stats, tGF, c("ms","seqgen")))

#for each f we will generate nreps alignments, and get all the stats for each
for (x in 1:length(fs)){
  
  print(x)
  
  for (y in 1:length(tGF)){

    receiver = paste("-es",tGF[y], "2", 1 - fs[x]) #f is 1 - proportion staying. (See ms documentation)
    donor = paste("-ej", tGF[y], "5 3")

    
    ms_opts = paste("-T",theta,pop_sizes,join1,join2,join3,receiver,donor,rec)
    
    sim_data_list <- mclapply(1:nreps, sim_D_f, nsam=nsam, P1=P1, P2=P2, P3=P3, P4=P4, P31=P31, P32=P32, ms_opts=ms_opts, len=len, scale=scale, subModel=subModel, stats=stats, mc.cores = 24)
    
    sim_data_array <- aaply(laply(sim_data_list, as.matrix),c(1,3),function(x) x) 
    
    all_data_32[x,,,y,] <- sim_data_array
            
    }
  
  }

all_data_32 <- ifelse(abs(all_data_32 == "Inf"), NA, all_data_32) # remove all values estimated calculated to be Inf or -Inf
  
all_data_32_D0 <- all_data_32
all_data_32_D0[,,"fG",,] <- ifelse(all_data_32[,,"D",,] >= 0, all_data_32[,,"fG",,], NA)
all_data_32_D0[,,"fhom",,] <- ifelse(all_data_32[,,"D",,] >= 0, all_data_32[,,"fhom",,], NA)
all_data_32_D0[,,"fd",,] <- ifelse(all_data_32[,,"D",,] >= 0, all_data_32[,,"fd",,], NA)

mean_data_32 <- apply(all_data_32,c(1,3,4,5),mean, na.rm = T) # average over all reps for each f
mean_data_32_D0 <- apply(all_data_32_D0,c(1,3,4,5),mean, na.rm = T) # average over all reps for each f

sd_data_32 <- apply(all_data_32,c(1,3,4,5),sd,na.rm=T) # get std err for all reps for each f
sd_data_32_D0 <- apply(all_data_32_D0,c(1,3,4,5),sd,na.rm=T) # get std err for all reps for each f


###gene flow 2-3
all_data_23 <- array(dim = c(length(fs), nreps, length(stats), length(tGF), length(c("ms","seqgen"))), dimnames = list(fs, 1:nreps, stats, tGF, c("ms","seqgen")))


for (x in 1:length(fs)){
  
  print(x)
  
  for (y in 1:length(tGF)){
    
    receiver = paste("-es", tGF[y], "3", 1 - fs[x]) #f is 1 - proportion staying. (See ms documentation)
    donor = paste("-ej", tGF[y], "5 2")
    
    
    ms_opts = paste("-T",theta,pop_sizes,join1,join2,join3,receiver,donor,rec)
    
    sim_data_list <- mclapply(1:nreps, sim_D_f, nsam=nsam, P1=P1, P2=P2, P3=P3, P4=P4, P31=P31, P32=P32, ms_opts=ms_opts, len=len, scale=scale, subModel=subModel, stats=stats, mc.cores = 24)
    
    sim_data_array <- aaply(laply(sim_data_list, as.matrix),c(1,3),function(x) x) 
    
    all_data_23[x,,,y,] <- sim_data_array
            
    }
  
  }

all_data_23 <- ifelse(abs(all_data_23 == "Inf"), NA, all_data_23) # remove all values estimated calculated to be Inf or -Inf
  
all_data_23_D0 <- all_data_23
all_data_23_D0[,,"fG",,] <- ifelse(all_data_23[,,"D",,] >= 0, all_data_23[,,"fG",,], NA)
all_data_23_D0[,,"fhom",,] <- ifelse(all_data_23[,,"D",,] >= 0, all_data_23[,,"fhom",,], NA)
all_data_23_D0[,,"fd",,] <- ifelse(all_data_23[,,"D",,] >= 0, all_data_23[,,"fd",,], NA)

mean_data_23 <- apply(all_data_23,c(1,3,4,5),mean, na.rm = T) # average over all reps for each f
mean_data_23_D0 <- apply(all_data_23_D0,c(1,3,4,5),mean, na.rm = T) # average over all reps for each f

sd_data_23 <- apply(all_data_23,c(1,3,4,5),sd,na.rm=T) # get std err for all reps for each f
sd_data_23_D0 <- apply(all_data_23_D0,c(1,3,4,5),sd,na.rm=T) # get std err for all reps for each f


######################################################################################
### plot data from seqgen simulations for Figure S1

jit = c(-0.01,0.01)
cols = c("red","blue")
#plot
png(file = paste("f_predictors", "_reps", nreps, "_l", len, "_r", rec_number, "_GF", paste(tGF, collapse = "_"), "_seqgen", "_scale", scale, ".png", sep = ""), width = 2000, height = 3500, res = 400)

par(mfrow = c(4,2), mar = c(3.5,3.5,2,1))

#D stat
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["3"]," to ","P"["2"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("D")), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_32[,"D",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_32[,"D",y,"seqgen"] - sd_data_32[,"D",y,"seqgen"], fs+jit[y], mean_data_32[,"D",y,"seqgen"] + sd_data_32[,"D",y,"seqgen"], col = cols[y])
  }

plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["2"]," to ","P"["3"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("D")), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_23[,"D",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_23[,"D",y,"seqgen"] - sd_data_23[,"D",y,"seqgen"], fs+jit[y], mean_data_23[,"D",y,"seqgen"] + sd_data_23[,"D",y,"seqgen"], col = cols[y])
  }

#f stat
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["3"]," to ","P"["2"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("f"["G"])), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_32_D0[,"fG",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_32_D0[,"fG",y,"seqgen"] - sd_data_32_D0[,"fG",y,"seqgen"], fs+jit[y], mean_data_32_D0[,"fG",y,"seqgen"] + sd_data_32_D0[,"fG",y,"seqgen"], col = cols[y])
  }
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["2"]," to ","P"["3"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("f"["G"])), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_23_D0[,"fG",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_23_D0[,"fG",y,"seqgen"] - sd_data_23_D0[,"fG",y,"seqgen"], fs+jit[y], mean_data_23_D0[,"fG",y,"seqgen"] + sd_data_23_D0[,"fG",y,"seqgen"], col = cols[y])
  }

#fhom stat
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["3"]," to ","P"["2"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("f"["hom"])), line = 2)
abline(0,1,lty=2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_32_D0[,"fhom",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_32_D0[,"fhom",y,"seqgen"] - sd_data_32_D0[,"fhom",y,"seqgen"], fs+jit[y], mean_data_32_D0[,"fhom",y,"seqgen"] + sd_data_32_D0[,"fhom",y,"seqgen"], col = cols[y])
  }
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["2"]," to ","P"["3"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("f"["hom"])), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_23_D0[,"fhom",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_23_D0[,"fhom",y,"seqgen"] - sd_data_23_D0[,"fhom",y,"seqgen"], fs+jit[y], mean_data_23_D0[,"fhom",y,"seqgen"] + sd_data_23_D0[,"fhom",y,"seqgen"], col = cols[y])
  }


#fd stat
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["3"]," to ","P"["2"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("f"["d"])), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_32_D0[,"fd",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_32_D0[,"fd",y,"seqgen"] - sd_data_32_D0[,"fd",y,"seqgen"], fs+jit[y], mean_data_32_D0[,"fd",y,"seqgen"] + sd_data_32_D0[,"fd",y,"seqgen"], col = cols[y])
  }
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, main = expression(paste("P"["2"]," to ","P"["3"])), xlab = "", ylab = "")
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(italic("f"["d"])), line = 2)
abline(0,1,lty=2)
for (y in 1:length(tGF)){
  points(fs+jit[y],mean_data_23_D0[,"fd",y,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19, col = cols[y])
  segments(fs+jit[y], mean_data_23_D0[,"fd",y,"seqgen"] - sd_data_23_D0[,"fd",y,"seqgen"], fs+jit[y], mean_data_23_D0[,"fd",y,"seqgen"] + sd_data_23_D0[,"fd",y,"seqgen"], col = cols[y])
  }

dev.off()




####################################################################################################
### simplified diagram for Figure 2


png(file = "Figure_2.png", width = 2000, height = 2000, res = 400)

par(mfrow = c(2,2), mar = c(2.25,2.25,0.25,0.25), ann = F)

#D stat
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, xaxt="n", yaxt="n", bty = "o")
abline(0,1,lty=2, lwd = 1.5, col = "gray")
points(fs,mean_data_32[,"D",1,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19)
segments(fs, mean_data_32[,"D",1,"seqgen"] - sd_data_32[,"D",1,"seqgen"], fs, mean_data_32[,"D",1,"seqgen"] + sd_data_32[,"D",1,"seqgen"], lwd = 1.5)
axis(2,at=c(0,.5,1))

plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, xaxt="n", yaxt="n", bty = "o")
abline(0,1,lty=2, lwd = 1.5, col = "gray")
points(fs,mean_data_23[,"D",1,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19)
segments(fs, mean_data_23[,"D",1,"seqgen"] - sd_data_23[,"D",1,"seqgen"], fs, mean_data_23[,"D",1,"seqgen"] + sd_data_23[,"D",1,"seqgen"], lwd = 1.5)


#fd
plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, xaxt="n", yaxt="n", bty = "o")
abline(0,1,lty=2, lwd = 1.5, col = "gray")
points(fs,mean_data_32[,"fd",1,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19)
segments(fs, mean_data_32[,"fd",1,"seqgen"] - sd_data_32[,"fd",1,"seqgen"], fs, mean_data_32[,"fd",1,"seqgen"] + sd_data_32[,"fd",1,"seqgen"], lwd = 1.5)
axis(1,at=c(0,.5,1))
axis(2,at=c(0,.5,1))

plot(0, xlim = c(0,1), ylim = c(0,1), cex = 0, xaxt="n", yaxt="n", bty = "o")
abline(0,1,lty=2, lwd = 1.5, col = "gray")
points(fs,mean_data_23[,"fd",1,"seqgen"], xlim = c(0,1), ylim = c(0,1), pch = 19)
segments(fs, mean_data_23[,"fd",1,"seqgen"] - sd_data_23[,"fd",1,"seqgen"], fs, mean_data_23[,"fd",1,"seqgen"] + sd_data_23[,"fd",1,"seqgen"], lwd = 1.5)
axis(1,at=c(0,.5,1))



dev.off()


