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

# Create the design file
design.df = data.frame(mutation = expDesign_df$Genotype,
					   row.names = rownames(expDesign_df)
					  )

# Create the DESeq2 matrix object and run DESeq
DESeq_obj = DESeqDataSetFromMatrix(countData = cts_m, colData = design.df, design = ~ mutation)
DESeq_obj = DESeq(DESeq_obj, fitType = "local")

# Get the results
mutationGenotypes = unique(expDesign_df$Genotype[which(expDesign_df$Genotype != "mrc1")])

# Create the result list
DESeq_resultsList = list()
for(i in 1:length(mutationGenotypes)){
	m = mutationGenotypes[i]
	DESeq_resultsList[[m]] = results(DESeq_obj, contrast = c("mutation", m, "mrc1"))
}

