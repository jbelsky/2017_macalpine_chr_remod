# Load the two files
brdu_cov_file = "/home/jab112/2017_macalpine_chr_remod/data/mrc1/DM409/DM409_brdu_coverage_2500bp_win.txt"
sim_file = "/home/jab112/2017_macalpine_chr_remod/data/mrc1/DM409/DM409_sim.txt"

brdu.df = read.table(brdu_cov_file, sep = "\t", header = T)
sim.df = read.table(brdu_cov_file, sep = "\t", header = T)

plot(density(log2(sim.df$Counts + 1)), type = "l")

# Obtain the cdf of the simulation
sim.cdf = ecdf(sim.v)

# Convert brdu coverage to quantiles based on simulation counts
quants.v = sim.cdf(brdu_cov.v)

# Get counts based on negative-binomial distribution
#counts = qnbinom(quants.v, size = 10, prob = 0.025)
