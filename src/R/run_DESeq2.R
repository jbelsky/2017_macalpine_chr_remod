library("DESeq2")

InitializeCtsMatrix <- function(ARG_expDesign_df, ARG_baseDir, ARG_fileSuffix){

	# Get the first sample coverage file
	expID = rownames(ARG_expDesign_df)[1]
	genotype = ARG_expDesign_df[expID, "Genotype"]
	coverageFile = paste(ARG_baseDir, "data", genotype, expID, paste0(expID, ARG_fileSuffix), sep = "/")

	# Read in the coverage file
	coverage_df = read.table(coverageFile, sep = "\t", header = T, row.names = 1)

	# Create the counts file
	ctsMatrix = matrix(0, nrow = nrow(coverage_df), ncol = nrow(ARG_expDesign_df))
	rownames(ctsMatrix) = rownames(coverage_df)
	colnames(ctsMatrix) = rownames(ARG_expDesign_df)

	# Return the ctsMatrix
	return(ctsMatrix)

}

GetIntegerCounts <- function(ARG_ctsMatrix, ARG_expDesign_df, ARG_expID, ARG_baseDir, ARG_fileSuffix){

	# Get the mutation genotype
	genotype = ARG_expDesign_df[ARG_expID, "Genotype"]

	# Get the coverage file
	coverageFile = paste(ARG_baseDir, "data", genotype, ARG_expID, paste0(ARG_expID, ARG_fileSuffix), sep = "/")

	# Read in the coverage files
	coverage_df = read.table(coverageFile, sep = "\t", header = T, row.names = 1)
	
	# Update the counts matrix
	ARG_ctsMatrix[, ARG_expID] = round(coverage_df$Norm_Counts)

	return(ARG_ctsMatrix)

}


# Set the file parameters
baseDir = system(command = "git rev-parse --show-toplevel", intern = T)
designFile = paste(baseDir, "GenomeFeatureFiles", "chr_remod_mutants.tsv", sep = "/")
fileSuffix = "_brdu_coverage_2500bp_win.txt"

# Load in the design dataframe
expDesign_df = read.table(designFile, sep = "\t", header = T, row.names = 1)

# Replace " " -> "_" in the "Genotype" column
expDesign_df$Genotype = gsub(pattern = " ", replacement = "_", x = expDesign_df$Genotype)

# Subset on the "double" mutants
expDesign_df = expDesign_df[which(expDesign_df$Mutant == "double"), ]

# Initialize the counts matrix
cts_m = InitializeCtsMatrix(expDesign_df, baseDir, fileSuffix)

# Iterate through each expID
for(r in 1:nrow(expDesign_df)){
	cts_m = GetIntegerCounts(cts_m, expDesign_df, rownames(expDesign_df)[r], baseDir, fileSuffix)
}

if(0){


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

}
