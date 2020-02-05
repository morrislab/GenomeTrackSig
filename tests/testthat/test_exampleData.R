# test_exampleData.R
# tests for running TrackSig end to end with example data

# ==== SETUP AND PREPARE =================================================

# loadable Hsapian genome
require(BSgenome.Hsapiens.UCSC.hg19)

# get the trajectory
vcfFile <- system.file("extdata", "Example.vcf", package = "TrackSig")
cnaFile <- system.file("extdata", "Example_cna.txt", package = "TrackSig")

detectedSigs <- detectActiveSignatures(vcfFile = vcfFile,
                       cnaFile = cnaFile,
                       purity = 1,
                       threshold = 0.035,
                       referenceSignatures = TrackSig:::alex_merged,
                       refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

traj <- TrackSig(sampleID = "example",
                 activeInSample = detectedSigs,
                 vcfFile = vcfFile,
                 cnaFile = cnaFile,
                 purity = 1,
                 referenceSignatures = TrackSig:::alex_merged,
                 scoreMethod = "SigFreq",
                 binSize = 100,
                 desiredMinSegLen = NULL,
                 refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

#plotTrajectory(traj)
#addPhiHist(traj, plotTrajectory(traj, linearX = T, anmac =F))
#
## ==== TEST ==============================================================
#
#test_that("a sample input prouces the expected output", {
#  expect_equal(T, T)
#})
#
#test_that("corrupt input generates errors",  {
#  expect_error(checkVcaf("a string", Hsapiens), "class(vcaf) not equal to \"data.frame\"", fixed = T)
#  expect_error(checkVcaf(vcaf, "a string"), "class(refGenome) not equal to \"BSgenome\"",  fixed = T)
#})
#

# ==== TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#


# [END]
