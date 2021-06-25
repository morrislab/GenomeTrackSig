library(ggplot2)

path <- "~/Desktop/CBSP2021/Cervix-SCC_pooled.csv"

detectedSigs <- detectActiveSigsCopy(path = path, binSize = 100)

traj <- TrackSigCopy(path = path, binSize = 100, activeInSample = detectedSigs, sampleID = "test")

plotTrajectory(traj, linearX = T) + labs(title = "Example trajectory with linear x-axis")

nonLinPlot <- plotTrajectory(traj, linearX = F) + labs(title = "Example trajectory with non-linear x-axis")
