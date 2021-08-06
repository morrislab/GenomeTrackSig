library(TrackSig)
library(ggplot2)

counts <- utils::read.csv(system.file(package = "TrackSig", "extdata/Example_counts.csv"), header=T, sep=',')

detectedSigs <- detectActiveSignatures(counts, threshold=0.15, binSize=300)

traj <- SpaceTrack(counts, activeInSample = detectedSigs, binSize=300)

plotSpaceTrajectory(traj, chr_level=F, cutoff=0) + labs("Example trajectory")

plotGenomicFeaturesTrajectory(counts, traj, 300, "Example trajectory with genomic features")
