# Get the "standard" sample

# Load the two files
sim_file = "/home/jab112/2017_macalpine_chr_remod/data/mrc1/DM409/DM409_sim.txt"
sim.df = read.table(sim_file, sep = "\t", header = T)

# Exclude chr 17 and rDNA on chr12
sim.df = sim.df[sim.df$Chromosome != 17,]
sim.df = sim.df[!(sim.df$Chromosome == 12 & sim.df$Position >= 450000 & sim.df$Position <= 490000),]

standard.v = sim.df$Counts

# Get the design file
design_file = "/home/jab112/2017_macalpine_chr_remod/GenomeFeatureFiles/chr_remod_mutants.tsv"
design.df = read.table(design_file, sep = "\t", header = T, stringsAsFactors = F)
data_dir = "/home/jab112/2017_macalpine_chr_remod/data"

for(dm_id in design.df[design.df$Mutant == "double","ID"]){
  
  print(dm_id)
  genotype = design.df[design.df$ID == dm_id, "Genotype"]
  genotype_str = gsub(" ", "_", genotype)
  
  sim_file = paste(data_dir, genotype_str, dm_id, paste0(dm_id, "_sim.txt"), sep = "/")
  sim.df = read.table(sim_file, sep = "\t", header = T)
  
  
  brdu_cov_file = paste(data_dir, genotype_str, dm_id, paste0(dm_id, "_brdu_coverage_2500bp_win.txt"), sep = "/")
  brdu.df = read.table(brdu_cov_file, sep = "\t", header = T)
  
  # Obtain the cdf of the simulation
  sim.cdf = ecdf(sim.df$Counts)
  
  # Convert brdu coverage to quantiles based on simulation counts
  quants.v = sim.cdf(brdu.df$Counts)
  
  # Get the normalized counts
  brdu.df[,"Norm_Counts"] = quantile(standard.v, probs = quants.v)
  
  write.table(brdu.df, file = brdu_cov_file, sep = "\t", col.names = T, row.names = F, quote = F)
  
}