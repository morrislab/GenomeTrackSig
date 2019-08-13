# test_callScripts.R
# tests for callScripts.R functions

# ==== SETUP AND PREPARE =================================================

# loadable Hsapian genome
#require(BSgenome.Hsapiens.UCSC.hg19)
#
#vcfFile <- system.file("extdata", "example.vcf", package = "TrackSig")
#vcaf <- getVcaf(vcfFile, cnaFile = NULL, purityFile = NULL, refGenome = Hsapiens)
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

