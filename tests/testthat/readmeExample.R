# readmeExample.R

library(GenomeTrackSig)
library(ggplot2)

counts <- utils::read.csv(system.file(package = "GenomeTrackSig", "extdata/Example_counts.csv"), header=T, sep=',')

detectedSigs <- detectActiveSignatures(counts, threshold=0.05, binSize=200)

traj <- GenomeTrackSig(counts, activeInSample = detectedSigs, binSize=200)

plotGenomeProfile(traj, chr_level=F, cutoff=0) + labs("Example profile")

plotGenomeContext(counts, traj, 200, "Example profile with additional genomic features")
