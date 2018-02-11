#!/usr/local/bin/R

library("TyphoonPlot")

# Load the BrdU csv
f_brdu = "/data/git/2017_macalpine_chr_remod/data/WT/DM278/DM278_WT_BrdU_OriDB_ACS798ori.csv"
f_brdu = "/data/git/2017_macalpine_chr_remod/data/mrc1/DM409/DM409_mrc1_BrdU_OriDB_ACS798ori.csv"
brdu.df = read.csv(file = f_brdu, header = T, sep = ",", stringsAsFactors = F, check.names = F)

# Convert to matrix
brdu.m = as.matrix(brdu.df[,-(1:4)])
rownames(brdu.m) = brdu.df[,1]

# Make the dens_dot_plot
png(file = "/home/jab112/public_html/20180211_plot.png",
    width = 14, height = 8, units = "in", res = 300
    )
par(mfrow = c(1,2))
par(mai = rep(1, 4))
DensDotPlot(brdu.m[order(rowMeans(brdu.m), decreasing = F)[1:250],])
par(mai = rep(1, 4))
plot(seq(-5000, 5000, 10), colMeans(brdu.m[order(rowMeans(brdu.m), decreasing = F)[1:250],]),
     type = "l", lwd = 2, col = "red"
)
dev.off()