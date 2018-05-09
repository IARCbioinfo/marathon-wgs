#!/usr/bin/env Rscript

library("Canopy")
# source("/home/pgm/Workspace/MPM/marathon/libs/custom_canopy.plottree.R")
source("/data/soudadel/MPM/canopy/scripts/custom_canopy.plottree.R")

##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

##
patient_id = args[1]
data_path  = args[2]

##
# patient_id = "5009"
# data_path  = "/home/pgm/Workspace/MPM/marathon/canopy/generated_clustering/5009/"

##
cat("####### ARGUMENTS #######\n")
cat(paste("patient_id: ", patient_id, "\n", sep=''))
cat(paste("output_path: ", data_path, "\n", sep=''))



##########################################
## Model selection with BIC
##########################################
bics = bicsn = c()
for (K in 2:9) {
  bic_tmp = try(read.table(file = paste(data_path, patient_id, ".K", K, '.BIC', sep="")), silent = TRUE)
  if (class(bic_tmp) != "try-error") {
    bics  = c(bics, bic_tmp)
    bicsn = c(bicsn, K)
  } else {
    cat(paste('BIC for K=', K, ' is not available.\n', sep=''))
  }
}
BIC = data.frame(as.numeric(bics))
bicsn = as.integer(bicsn)
BIC = cbind(bicsn, BIC)
rownames(BIC) = bicsn
colnames(BIC) = c("K", "BIC")
optimalBIC = BIC[which.max(BIC$BIC), ]


##########################################
## Plot K selection
##########################################
svg(filename=paste(data_path, patient_id, ".optimal_K.svg", sep=''))
plot(BIC, xlab="K (nb subclones)", xaxt="n")
lines(BIC)
abline(v = optimalBIC$K, col="red")
axis(side=1, at=bicsn, labels=bicsn)
legend("topleft", c(paste("optimal K = ", optimalBIC$K, sep='')), col=c(rgb(1,0,0,1)), lty=1:2, cex=0.8, box.lty=0)
dev.off() # close the device


##########################################
## Posterior evaluation of sampled trees, output, and plotting
##########################################
# Canopy then runs a posterior evaluation of all sampled trees by MCMC. 
# If modes of posterior probabilities (second column of config.summary) aren’t obvious, 
# check if the algorithm has converged (and run sampling longer if not).
load(file = paste(data_path, patient_id, ".K", optimalBIC$K, ".postmcmc_image.rda", sep=''))
post = canopy.post(sampchain = sampchain, 
                   projectname = data_path, 
                   K = 1,
                   numchain = numchain, 
                   burnin = 90, 
                   thin = 5, 
                   optK = 1, # because of parallelization, it's always [[1]]
                   C = C, 
                   post.config.cutoff = 0.05)

samptreethin     = post[[1]] # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]] # likelihoods of trees in samptree
config           = post[[3]] # configuration for each posterior tree
config.summary   = post[[4]] # configuration summary
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space


##########################################
## Tree
##########################################
# One can then use Canopy to output and plot the most likely tree (i.e.,tree with 
# the highest posterior likelihood). Mutations, clonal frequencies, and tree topology, 
# etc., of the tree are obtained from the posterior distributions of subtree space 
# with trees having the same configuration. In our MDA231 example, the most likely 
# tree is the tree with configuration 3. Note: A separate txt file can be generated 
# (with txt=TRUE and txt.name=’*.txt’) if the figure legend of mutational profiles 
# (texts below the phylogenetic tree) in the plot is too long to be fitted entirely.

config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood!\n')
# plot the most likely tree in the posterior tree space
output.tree = canopy.output(post, config.i, C)
custom_canopy.plottree(output.tree, pdf=TRUE, pdf.name = paste(data_path, patient_id, '_TREE_most_likely.pdf', sep=''), txt = TRUE, txt.name = paste(data_path, patient_id, '_TREE_most_likely.txt', sep=''))
