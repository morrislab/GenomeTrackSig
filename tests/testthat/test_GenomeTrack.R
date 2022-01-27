# library(GenomeTrackSig)
# library(ggplot2)
#
# counts <- utils::read.csv(system.file(package = "GenomeTrackSig", "extdata/Example_counts.csv"), header=T, sep=',')
#
# detectedSigs <- detectActiveSignatures(counts, threshold=0.15, binSize=300)
#
# traj <- GenomeTrackSig(counts, activeInSample = detectedSigs, binSize=300)
#
# plotGenomeProfile(traj, chr_level=F, cutoff=0) + labs("Example trajectory")
#
# plotGenomeContext(counts, traj, 300, "Example trajectory with genomic features")
