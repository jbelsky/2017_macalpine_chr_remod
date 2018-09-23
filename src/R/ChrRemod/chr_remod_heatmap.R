library("TyphoonPlot")

d_data = "/data/git/2017_macalpine_chr_remod/data"

chr_remod = c("mrc1", "mrc1_chd1_asf1", "mrc1_chd1_swr1", "mrc1_isw1_asf1",
              "mrc1_isw1_chd1", "mrc1_isw1_isw2", "mrc1_isw2_asf1", "mrc1_isw2_chd1",
              "mrc1_isw2_swr1"
             )

brdu.m = matrix(0, nrow = 798, ncol = length(chr_remod))

mat_idx = 1
for (i in chr_remod){
  
  # Get the tsv files
  f_rpkm = list.files(paste0(d_data, "/", i), pattern = ".tsv", recursive = T)
  a.df = read.table(paste0(d_data, "/", i, "/", f_rpkm[1]), header = T)
  b.df = read.table(paste0(d_data, "/", i, "/", f_rpkm[2]), header = T)
  
  # Average the rpkm into the brdu.m
  brdu.m[,mat_idx] = (a.df$rpkm + b.df$rpkm) / 2
  mat_idx = mat_idx + 1
  
  # Compare the replicates
  png(file = paste0("/home/jab112/public_html/2018_09_17/replicate_comparison/", i, ".png"),
      width = 7, height = 7, units = "in", res = 300
  )
  
  plot(a.df$rpkm, b.df$rpkm, cex = 0.5, pch = 19, 
       xlim = c(0, 500), ylim = c(0, 500),
       xlab = substring(f_rpkm[1], 1, 5),
       ylab = substring(f_rpkm[2], 1, 5),
       main = i
  )
  abline(0, 1, col = "red", lty = 2)
  
  dev.off()
  
}

for(i in 2:ncol(brdu.m)){

  cur_chr_remod = chr_remod[i]
  
  # Compare to mrc1
  png(file = paste0("/home/jab112/public_html/2018_09_17/comparison_to_mrc1/mrc1_to_", cur_chr_remod, ".png"),
      width = 7, height = 7, units = "in", res = 300
  )
  
  plot(brdu.m[,1], brdu.m[,i], cex = 0.5, pch = 19, 
       xlim = c(0, 500), ylim = c(0, 500),
       xlab = "mrc1",
       ylab = cur_chr_remod,
       main = paste0("mrc1 vs. ", cur_chr_remod)
  )
  abline(0, 1, col = "red", lty = 2)
  
  dev.off()

}

#early_idx = which(acs.df$activation_time == "early")
#brdu.m = brdu.m[early_idx,]
if(0){


idx = order(rowMeans(brdu.m), decreasing = F)
brdu.m = brdu.m[idx,]

DensDotPlot(brdu.m, z_max = 500)

mat.m = matrix(0, nrow = nrow(brdu.m), ncol = ncol(brdu.m))

for(i in 2:ncol(brdu.m)){
  
  # Get the log2 diff
  log2_diff.v = log2(brdu.m[,1] + 0.1) - log2(brdu.m[,i] + 0.1)
  loss.v = which(log2_diff.v <= -1)
  gain.v = which(log2_diff.v >= 1)
  mat.m[loss.v, i] = -1
  mat.m[gain.v, i] = 1
  
  mat.m[,i] = log2_diff.v
  
  
}

DensDotPlot(mat.m, z_min = -2, z_max = 2, low_col = "red", med_col = "white", high_col = "darkgreen")

}