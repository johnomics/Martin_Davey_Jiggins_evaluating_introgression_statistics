#!/usr/bin/env Rscript

# Figure_4.R
# Plot for simulations of D, f and dxy

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# August 2014


big_table <- read.delim("model_files_win10000_s0.01_l5000_r50.alternate_models.dxy.summary.sg.tsv", sep = "\t")

null_table <- read.delim("model_files_win10000_s0.01_l5000_r50.null_models.dxy.summary.sg.tsv", sep = "\t")

big_table_5 <- read.delim("model_files_win10000_s0.01_l5000_r5.alternate_models.dxy.summary.sg.tsv", sep = "\t")

null_table_5 <- read.delim("model_files_win10000_s0.01_l5000_r5.null_models.dxy.summary.sg.tsv", sep = "\t")


png("Figure_4.png", width = 3000, height = 1500, res = 400)
par(mfrow = c(2,4), mar = c(3,3.5,1,0.5))

pos = c(1,2,4,5,7,8)
pos_null = c(4,5,7,8)


bt123 = 0.8
bt21 = 0.6
at123 = 0.6
at23 = 0.4

plot_table <- subset(big_table, Background_t123 == bt123 & Background_t21 == bt21 & Alternate_t123 == at123 & Alternate_t23 == at23 & Variable == "P2P3_dxy")

plot(pos,plot_table$Mean[c(1,2,3,4,9,10)], ylim = c(0,0.033), col = c("#ad5c08", "#00ba38", rep(c( "gray", "black"),4)), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos,plot_table$Mean[c(1,2,3,4,9,10)]+plot_table$SD[c(1,2,3,4,9,10)], pos, plot_table$Mean[c(1,2,3,4,9,10)]-plot_table$SD[c(1,2,3,4,9,10)], col = c("#ad5c08", "#00ba38", rep(c( "gray", "black"),4)))
axis(1, at = c(1.5,4.5, 7.5), labels = c("Sim.", expression(italic("D")), expression(italic("f"["d"]))))
mtext(side = 2, text = expression(italic("d"["XY"])), line = 2)
text(c(1.5,4.5,7.5), rep(0.01,3), ifelse(plot_table$w.p.background.higher[c(1,3,9)] < 0.001, "*", ""), cex = 2)

bt123 = 0.8
bt21 = 0.6
at123 = 0.8
at23 = 0.4

plot_table <- subset(big_table, Background_t123 == bt123 & Background_t21 == bt21 & Alternate_t123 == at123 & Alternate_t23 == at23 & Variable == "P2P3_dxy")

plot(pos,plot_table$Mean[c(1,2,3,4,9,10)], ylim = c(0,0.033), col = c("#ad5c08", "#619cff", rep(c( "gray", "black"),4)), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos,plot_table$Mean[c(1,2,3,4,9,10)]+plot_table$SD[c(1,2,3,4,9,10)],pos, plot_table$Mean[c(1,2,3,4,9,10)]-plot_table$SD[c(1,2,3,4,9,10)], col = c("#ad5c08", "#619cff", rep(c( "gray", "black"),4)))
axis(1, at = c(1.5,4.5, 7.5), labels = c("Sim.", expression(italic("D")), expression(italic("f"["d"]))))
text(c(1.5,4.5,7.5), rep(0.01,3), ifelse(plot_table$w.p.background.higher[c(1,3,9)] < 0.001, "*", ""), cex = 2)



bt123 = 0.8
bt21 = 0.6
at123 = 1
at23 = 0.8

plot_table <- subset(big_table, Background_t123 == bt123 & Background_t21 == bt21 & Alternate_t123 == at123 & Alternate_t23 == at23 & Variable == "P2P3_dxy")

plot(pos,plot_table$Mean[c(1,2,3,4,9,10)], ylim = c(0,0.033), col = c("#ad5c08", "#f8766d", rep(c( "gray", "black"),4)), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos,plot_table$Mean[c(1,2,3,4,9,10)]+plot_table$SD[c(1,2,3,4,9,10)], pos, plot_table$Mean[c(1,2,3,4,9,10)]-plot_table$SD[c(1,2,3,4,9,10)], col = c("#ad5c08", "#f8766d", rep(c( "gray", "black"),4)))
axis(1, at = c(1.5,4.5, 7.5), labels = c("Sim.", expression(italic("D")), expression(italic("f"["d"]))))
text(c(1.5,4.5,7.5), rep(0.01,3), ifelse(plot_table$w.p.background.higher[c(1,3,9)] < 0.001, "*", ""), cex = 2)




bt123 = 0.8
bt21 = 0.6

plot_table <- subset(null_table, Background_t123 == bt123 & Background_t21 == bt21 & Variable == "P2P3_dxy")

plot(pos_null,plot_table$Mean[c(1,2,7,8)], ylim = c(0,0.033), col = rep(c( "gray", "black"),5), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos_null,plot_table$Mean[c(1,2,7,8)]+plot_table$SD[c(1,2,7,8)], pos_null, plot_table$Mean[c(1,2,7,8)]-plot_table$SD[c(1,2,7,8)], col = rep(c( "gray", "black"),5))
axis(1, at = c(4.5, 7.5), labels = c(expression(italic("D")), expression(italic("f"["d"]))))
text(c(4.5,7.5), rep(0.01,2), ifelse(plot_table$w.p.background.higher[c(1,7)] < 0.001, "*", ""), cex = 2)


### low rec


bt123 = 0.8
bt21 = 0.6
at123 = 0.6
at23 = 0.4

plot_table <- subset(big_table_5, Background_t123 == bt123 & Background_t21 == bt21 & Alternate_t123 == at123 & Alternate_t23 == at23 & Variable == "P2P3_dxy")

plot(pos,plot_table$Mean[c(1,2,3,4,9,10)], ylim = c(0,0.033), col = c("#ad5c08", "#00ba38", rep(c( "gray", "black"),4)), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos,plot_table$Mean[c(1,2,3,4,9,10)]+plot_table$SD[c(1,2,3,4,9,10)], pos, plot_table$Mean[c(1,2,3,4,9,10)]-plot_table$SD[c(1,2,3,4,9,10)], col = c("#ad5c08", "#00ba38", rep(c( "gray", "black"),4)))
axis(1, at = c(1.5,4.5, 7.5), labels = c("Sim.", expression(italic("D")), expression(italic("f"["d"]))))
mtext(side = 2, text = expression(italic("d"["XY"])), line = 2)
text(c(1.5,4.5,7.5), rep(0.01,3), ifelse(plot_table$w.p.background.higher[c(1,3,9)] < 0.001, "*", ""), cex = 2)


bt123 = 0.8
bt21 = 0.6
at123 = 0.8
at23 = 0.4

plot_table <- subset(big_table_5, Background_t123 == bt123 & Background_t21 == bt21 & Alternate_t123 == at123 & Alternate_t23 == at23 & Variable == "P2P3_dxy")

plot(pos,plot_table$Mean[c(1,2,3,4,9,10)], ylim = c(0,0.033), col = c("#ad5c08", "#619cff", rep(c( "gray", "black"),4)), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos,plot_table$Mean[c(1,2,3,4,9,10)]+plot_table$SD[c(1,2,3,4,9,10)],pos, plot_table$Mean[c(1,2,3,4,9,10)]-plot_table$SD[c(1,2,3,4,9,10)], col = c("#ad5c08", "#619cff", rep(c( "gray", "black"),4)))
axis(1, at = c(1.5,4.5, 7.5), labels = c("Sim.", expression(italic("D")), expression(italic("f"["d"]))))
text(c(1.5,4.5,7.5), rep(0.01,3), ifelse(plot_table$w.p.background.higher[c(1,3,9)] < 0.001, "*", ""), cex = 2)



bt123 = 0.8
bt21 = 0.6
at123 = 1
at23 = 0.8

plot_table <- subset(big_table_5, Background_t123 == bt123 & Background_t21 == bt21 & Alternate_t123 == at123 & Alternate_t23 == at23 & Variable == "P2P3_dxy")

plot(pos,plot_table$Mean[c(1,2,3,4,9,10)], ylim = c(0,0.033), col = c("#ad5c08", "#f8766d", rep(c( "gray", "black"),4)), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos,plot_table$Mean[c(1,2,3,4,9,10)]+plot_table$SD[c(1,2,3,4,9,10)], pos, plot_table$Mean[c(1,2,3,4,9,10)]-plot_table$SD[c(1,2,3,4,9,10)], col = c("#ad5c08", "#f8766d", rep(c( "gray", "black"),4)))
axis(1, at = c(1.5,4.5, 7.5), labels = c("Sim.", expression(italic("D")), expression(italic("f"["d"]))))
text(c(1.5,4.5,7.5), rep(0.01,3), ifelse(plot_table$w.p.background.higher[c(1,3,9)] < 0.001, "*", ""), cex = 2)




bt123 = 0.8
bt21 = 0.6

plot_table <- subset(null_table_5, Background_t123 == bt123 & Background_t21 == bt21 & Variable == "P2P3_dxy")

plot(pos_null,plot_table$Mean[c(1,2,7,8)], ylim = c(0,0.033), col = rep(c( "gray", "black"),5), pch = 19, xlim = c(0,9),yaxp = c(0,0.03,3), xaxt = "n", ylab = "")
segments(pos_null,plot_table$Mean[c(1,2,7,8)]+plot_table$SD[c(1,2,7,8)], pos_null, plot_table$Mean[c(1,2,7,8)]-plot_table$SD[c(1,2,7,8)], col = rep(c( "gray", "black"),5))
axis(1, at = c(4.5, 7.5), labels = c(expression(italic("D")), expression(italic("f"["d"]))))
text(c(4.5,7.5), rep(0.01,2), ifelse(plot_table$w.p.background.higher[c(1,7)] < 0.001, "*", ""), cex = 2)



dev.off()




