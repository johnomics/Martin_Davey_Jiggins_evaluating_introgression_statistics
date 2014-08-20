#!/usr/bin/env Rscript

# Figure_1.R
# Use the derivation for expected D (Durand et al. 2011) to plot the effects of Ne.

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# August 2014


####################################################################################
####### D as a function of f

png(file = "Figure_1B.png", width = 2000, height = 2000, res = 400)
par(mar=c(4,4,1,1), ann = F)


N = 500000

tp2 <- 1000000

tp3 <- 2000000

tgf = 0


f = seq(0,1,0.01)

D_exp <- (3*f*(tp3 - tgf)) / (3*f*(tp3 - tgf) + 4*N*(1-f)*(1-1/(2*N))^(tp3-tp2) + 4*N*f*(1-1/(2*N))^(tp3 - tgf))

plot(D_exp ~ f, type = "l", ylim = c(0,1), lty = 1)
mtext(side = 1, text = expression(italic("f")), line = 2)
mtext(side = 2, text = expression(paste(italic("E"),"[",italic("D"),"]")), line = 2)




N = 1000000

tp2 <- 1000000

tp3 <- 2000000

tgf = 0


f = seq(0,1,0.01)

D_exp <- (3*f*(tp3 - tgf)) / (3*f*(tp3 - tgf) + 4*N*(1-f)*(1-1/(2*N))^(tp3-tp2) + 4*N*f*(1-1/(2*N))^(tp3 - tgf))

lines(D_exp ~ f, type = "l", ylim = c(0,1), lty = 2)



N = 2000000

tp2 <- 1000000

tp3 <- 2000000

tgf = 0


f = seq(0,1,0.01)

D_exp <- (3*f*(tp3 - tgf)) / (3*f*(tp3 - tgf) + 4*N*(1-f)*(1-1/(2*N))^(tp3-tp2) + 4*N*f*(1-1/(2*N))^(tp3 - tgf))

lines(D_exp ~ f, type = "l", ylim = c(0,1), lty = 3)



legend(0.7, 0.3, title = expression(paste("N"["e"]," (million)")),legend = c(0.5, 1, 2), lty = c(1, 2, 3), bty = "n") 

dev.off()




#####################################################################################################
### D as a function of Ne

png(file = "Figure_1C.png", width = 2000, height = 2000, res = 400)

par(mar=c(4,4,1,1), ann = F)

N = seq(10000,2000000,10000)

tp2 <- 500000

tp3 <- 2000000

tgf = 0


f = 0.1

D_exp <- (3*f*(tp3 - tgf)) / (3*f*(tp3 - tgf) + 4*N*(1-f)*(1-1/(2*N))^(tp3-tp2) + 4*N*f*(1-1/(2*N))^(tp3 - tgf))

plot(D_exp ~ N, type = "l", ylim = c(0,1))
mtext(side = 2, text = expression(paste(italic("E"),"[",italic("D"),"]")), line = 2)
mtext(side = 1, text = expression("N"["e"]), line = 2)




tp2 <- 1000000

tp3 <- 2000000

tgf = 0


f = 0.1

D_exp <- (3*f*(tp3 - tgf)) / (3*f*(tp3 - tgf) + 4*N*(1-f)*(1-1/(2*N))^(tp3-tp2) + 4*N*f*(1-1/(2*N))^(tp3 - tgf))

lines(D_exp ~ N, type = "l", lty = 2)




tp2 <- 1500000

tp3 <- 2000000

tgf = 0


f = 0.1

D_exp <- (3*f*(tp3 - tgf)) / (3*f*(tp3 - tgf) + 4*N*(1-f)*(1-1/(2*N))^(tp3-tp2) + 4*N*f*(1-1/(2*N))^(tp3 - tgf))

lines(D_exp ~ N, type = "l", lty = 3)


legend(750000, 1, title = expression(paste(italic("t")["12"], " (million generations ago)")), legend = c(0.5,1,1.5), lty = c(1, 2, 3), bty = "n") 


dev.off()
