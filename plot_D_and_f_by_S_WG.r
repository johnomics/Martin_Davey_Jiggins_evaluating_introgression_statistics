#!/usr/bin/env Rscript

# plot_D_and_f_by_S.r
# Plots for whole genome Heliconius data - Figure 5

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# November-December 2013



############## Functions ###############

#function to rearrange scaffolds by chromosome and add a chromosome column
as.chromosomes <- function(table,agp,chromNames,gap = 0) {
  WG_table <- table[0,]
  chromosome <- vector()
  genome_pos <- vector()
  chrom_pos <- vector()
  i = 1
  endLast = 0
  for (name in chromNames){
    agpCurrent <- subset(agp, chromosome == name)
    starts = as.numeric(agpCurrent$start)
    ends = as.numeric(agpCurrent$end)
    chrom_table <- table[0,]
    for (x in 1:length(agpCurrent$scaffold)) {
      scaf_table <- subset(table, scaffold == agpCurrent$scaffold[x])
      chrom_table <- rbind(chrom_table,scaf_table)
      for (y in as.integer(scaf_table$position)) {
        if (agpCurrent$ori[x] == "+"){
          chrom_pos[i] <- y + starts[x]
          genome_pos[i] <- y + starts[x] + endLast + gap
          }
        else {
          chrom_pos[i] <- ends[x] - y
          genome_pos[i] <- ends[x] - y  + endLast + gap
          }
        chromosome[i] <- name
        i <- i + 1  
        }
      }
    WG_table <- rbind(WG_table,chrom_table)
    endLast = endLast + gap + ends[length(ends)]
    }
  WG_table <- cbind(chromosome,WG_table,chrom_pos,genome_pos)
  WG_table <- WG_table[with(WG_table,order(genome_pos)),]
  return(WG_table)
  }


########### input data ################

#agp file carrys chromosome information for mapped scaffolds
agp <- read.delim("Hmel1-1_hox_RAD_matepair_chromosomes_Zupdated.agp", as.is = T, sep = "\t", header = F)

names(agp) = c("chromosome", "start", "end", "order", "DN", "scaffold", "one", "length", "ori", "rec")
agp <- subset(agp, DN == "D")


chromNames <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chrZ")


#input data

ABBA_table <- read.csv("set31.Zupdated.autoscafs.union.SPiDxyDf.w5m1s5.csv", as.is = T)
ABBA_table_Z <- read.csv("set31.Zupdated.chrZscafs.union.SPiDxyDf.w5m1s5.csv", as.is = T)
ABBA_table <- rbind(ABBA_table, ABBA_table_Z)


#filter for windows with enough sites above the minimum data cutoff
ABBA_table <- subset(ABBA_table, sitesOverMinExD >= 1000)

#rearrange scaffolds by chromosome (this is only necessary for the third, chromosomal, plot)
WG_table <- as.chromosomes(ABBA_table,agp,chromNames)
WG_table <- subset(WG_table, chromosome %in% chromNames)

#scaffolds containing wing patterning loci
ybsub <- subset(ABBA_table, scaffold == "HE667780" & position >= 600000 & position <= 1000000)
bdsub <- subset(ABBA_table, scaffold == "HE670865" & position >= 300000 & position <= 500000)

###### plotting D and f against number of segregating sites


pdf(file = "D_and_f_against_s.pdf", width = 3.3, height = 9.5)

par(mfrow = c(3,1), mar = c(3,3,1,1))

plot(ABBA_table$S/ABBA_table$sitesOverMinExD, ABBA_table$D, cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.05), bty = "n", xlab = "", ylab = "", xlim = c(0.05,0.35) )
points(ybsub$S/ybsub$sitesOverMinExD, ybsub$D, cex = 1, pch = 1, col = rgb(1,0.9,0,1))
points(bdsub$S/bdsub$sitesOverMinExD, bdsub$D, cex = 1, pch = 1, col = "red")
segments(mean(ABBA_table$S/ABBA_table$sitesOverMinExD, na.rm = T), -1, mean(ABBA_table$S/ABBA_table$sitesOverMinExD, na.rm = T), 1, lty = 3)
abline(quantile(ABBA_table$D,0.9,na.rm=T), 0, lty = 3)
abline(quantile(ABBA_table$D,0.1,na.rm=T), 0, lty = 3)


plot(ABBA_table$S[ABBA_table$D >= 0]/ABBA_table$sitesOverMinExD[ABBA_table$D >= 0], ABBA_table$mf[ABBA_table$D >= 0], cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.05), ylim = c(0,1), bty = "n", xlab = "", ylab = "", xlim = c(0.05,0.35) )
points(ybsub$S[ybsub$D >= 0]/ybsub$sitesOverMinExD[ybsub$D >= 0], ybsub$mf[ybsub$D >= 0], cex = 1, pch = 1, col = rgb(1,0.9,0,1))
points(bdsub$S[bdsub$D >= 0]/bdsub$sitesOverMinExD[bdsub$D >= 0], bdsub$mf[bdsub$D >= 0], cex = 1, pch = 1, col = "red")
segments(mean(ABBA_table$S[ABBA_table$D >= 0]/ABBA_table$sitesOverMinExD[ABBA_table$D >= 0], na.rm = T), 0, mean(ABBA_table$S[ABBA_table$D >= 0]/ABBA_table$sitesOverMinExD[ABBA_table$D >= 0], na.rm = T), 1, lty = 3)
abline(quantile(ABBA_table$mf[ABBA_table$D >= 0],0.9,na.rm=T), 0, lty = 3)


### variance by number of segregating sites - plotted by chromosome
library(RColorBrewer)
cols <- c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"))


#means and variances
vD <- numeric(length <- length(chromNames))
vmf <- numeric(length <- length(chromNames))
meanS <- numeric(length <- length(chromNames))

for (x in 1:length(chromNames)){
  vD[x] <- var(WG_table$D[WG_table$chromosome == chromNames[x]], na.rm = T)
  vmf[x] <- var(WG_table$mf[WG_table$chromosome == chromNames[x] & WG_table$D >= 0], na.rm = T)
  meanS[x] <- mean(WG_table$S[WG_table$chromosome == chromNames[x]] / WG_table$sitesOverMinExD[WG_table$chromosome == chromNames[x]], na.rm = T)
  }


plot(vD ~ meanS, ylim = c(0,0.16), xlim = c(0.12,0.18), pch = 19, cex = 1.5, col = cols, xlab = "", ylab = "", bty = "n")

points(vmf ~ meanS, pch = 15, cex = 1.5, col = cols)

legend(0.16, 0.15,legend = c(1:20,"Z"), pch = 19, col = cols, ncol = 3, bty = "n", cex = 1, pt.cex = 1.2)

dev.off()
