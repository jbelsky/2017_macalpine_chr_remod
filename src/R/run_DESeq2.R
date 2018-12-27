library("DESeq2")

# Read in the coverage files
dm409.df = read.table("/home/jab112/2017_macalpine_chr_remod/data/mrc1/DM409/DM409_brdu_coverage_2500bp_win.txt", sep = "\t", header = T, row.names = 1)
dm410.df = read.table("/home/jab112/2017_macalpine_chr_remod/data/mrc1/DM410/DM410_brdu_coverage_2500bp_win.txt", sep = "\t", header = T, row.names = 1)
dm391.df = read.table("/home/jab112/2017_macalpine_chr_remod/data/mrc1_isw1_isw2/DM391/DM391_brdu_coverage_2500bp_win.txt", sep = "\t", header = T, row.names = 1)
dm414.df = read.table("/home/jab112/2017_macalpine_chr_remod/data/mrc1_isw1_isw2/DM414/DM414_brdu_coverage_2500bp_win.txt", sep = "\t", header = T, row.names = 1)

# Create the counts file
cts.m = matrix(0, nrow = nrow(dm409.df), ncol = 4)
rownames(cts.m) = rownames(dm409.df)
cts.m[,1] = round(dm409.df$Norm_Counts)
cts.m[,2] = round(dm410.df$Norm_Counts)
cts.m[,3] = round(dm391.df$Norm_Counts)
cts.m[,4] = round(dm414.df$Norm_Counts)
colnames(cts.m) = c("ctrl1", "ctrl2", "mutant1", "mutant2")

# Create the design file
design.df = data.frame(mutation = c("ctrl", "ctrl", "mutant", "mutant"),
					   row.names = colnames(cts.m)
					  )

# Create the DESeq2 matrix object and run DESeq
dds = DESeqDataSetFromMatrix(countData = cts.m, colData = design.df, design = ~ mutation)
dds = DESeq(dds, fitType = "local")

# Get the results
res = results(dds)
result.df = as.data.frame(res)
