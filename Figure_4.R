#!/usr/bin/env Rscript

# Figure_4.R
# Generate plots across BD and Yb bac walk regions for Figure 4 as a PDF file

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013


#input data

yb <- read.csv("Heliconius_HmB_windows.csv")
yb <- subset(yb, is.na(D) == F & is.na(mf) == F & sitesOverMinExD >= 1000)


bd <- read.csv("Heliconius_HmYb_windows.csv")
bd <- subset(bd, is.na(D) == F & is.na(mf) == F & sitesOverMinExD >= 1000)

pdf("Figure_4.pdf", width = 12, height = 7)


par(mfrow = c(4,2), mar = c(2.5,2,0.5,1))

#counts of ABBA and BABA per 100 sites
plot(0,0,cex = 0, ylim = c(0,1.2), bty = "n", xaxt = "n", xlab = "", xlim = c(0,1200000), ylab = "count per 100 sites")
polygon(c(min(yb$position),yb$position,max(yb$position)),c(0,yb$ABBA*100/yb$sites,0), border = NA, col = rgb(255,102,164,150, maxColorValue = 255))
polygon(c(min(yb$position),yb$position,max(yb$position)),c(0,yb$BABA*100/yb$sites,0), border = NA, col = rgb(60,192,176,150, maxColorValue = 255))

legend(10000, 1.2, legend = c("ABBA", "BABA"), fill = c(rgb(255,102,164,150, maxColorValue = 255),rgb(60,192,176,150, maxColorValue = 255)), bty = "n", border = NA, cex = 0.9)

plot(0,0,cex = 0, ylim = c(0,1.2), bty = "n", xaxt = "n", xlab = "", xlim = c(0,1200000), ylab = "count per 100 sites")
polygon(c(min(bd$position),bd$position,max(bd$position)),c(0,bd$ABBA*100/bd$sites,0), border = NA, col = rgb(255,102,164,150, maxColorValue = 255))
polygon(c(min(bd$position),bd$position,max(bd$position)),c(0,bd$BABA*100/bd$sites,0), border = NA, col = rgb(60,192,176,150, maxColorValue = 255))


# D stat
plot(0,0,cex = 0, ylim = c(-1,1), bty = "n", xlab = "", xlim = c(0,1200000), ylab = "f", xaxt = "n")
polygon(c(min(yb$position),yb$position,max(yb$position)),c(0,yb$D,0), border = NA, col = rgb(0,0,0,0.5))

plot(0,0,cex = 0, ylim = c(-1,1), bty = "n", xlab = "", xlim = c(0,1200000), ylab = "f", xaxt = "n")
polygon(c(min(bd$position),bd$position,max(bd$position)),c(0,bd$D,0), border = NA, col = rgb(0,0,0,0.5))



# f stat
plot(0,0,cex = 0, ylim = c(0,1), bty = "n", xlab = "", xlim = c(0,1200000), ylab = "f", xaxt = "n")
polygon(c(min(yb$position[yb$D >= 0]),yb$position[yb$D >= 0],max(yb$position[yb$D >= 0])),c(0,yb$mf[yb$D >= 0],0), border = NA, col = rgb(0,0,0,0.5))

plot(0,0,cex = 0, ylim = c(0,1), bty = "n", xlab = "", xlim = c(0,1200000), ylab = "f", xaxt = "n")
polygon(c(min(bd$position[bd$D >= 0]),bd$position[bd$D >= 0],max(bd$position[bd$D >= 0])),c(0,bd$mf[bd$D >= 0],0), border = NA, col = rgb(0,0,0,0.5))


# dXY

plot(0,0,cex = 0, ylim = c(0,0.06), bty = "n", xlab = "", xlim = c(0,1200000), ylab = "f", xaxt = "n")
polygon(c(min(yb$position),yb$position[is.na(yb$P2P3_dxy) == F],max(yb$position)),c(0,yb$P2P3_dxy[is.na(yb$P2P3_dxy) == F],0), border = NA, col = rgb(0,0,0,0.5))
axis(1,at = seq(0,1200000,200000), labels = seq(0, 1200, 200))

plot(0,0,cex = 0, ylim = c(0,0.06), bty = "n", xlab = "", xlim = c(0,1200000), ylab = "f", xaxt = "n")
polygon(c(min(bd$position),bd$position[is.na(bd$P2P3_dxy) == F],max(bd$position)),c(0,bd$P2P3_dxy[is.na(bd$P2P3_dxy) == F],0), border = NA, col = rgb(0,0,0,0.5))
axis(1,at = seq(0,800000,200000), labels = seq(0,800,200))


dev.off()

