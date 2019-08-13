# call_scripts.R


#' @export
vcfToCounts_simulation <- function(vcfFile, mutTypesFile, cnaFile = NULL, purityFile = NULL,
                        context = generateContext(c("CG", "TA")), refGenome = BSgenome.Hsapiens.UCSC.hg19, binSize = 100,
                        saveIntermediate = F, intermediateFile = NULL) {

  # load CNA and purity dataframe (not loaded with VCF for parallelization memory saving)
  # could be done as a single annotation load.... one function to load each file
  # loads the following - all shared between all VCF's, all optional (but not necessarily independent)
  # cna, purity, tumortypes, signatures (alex, cosmic), trinucleotide, sigactivities
  # tumortype_file = "", signature_file = "", trinucleotide_file = "", active_signatures_file = ""

  # input checking and path expansion


  # input checking and path expansion
  vcfFile <- path.expand(vcfFile)
  stopifnot(file.exists(vcfFile))

  if (!is.null(cnaFile)){
    cnaFile <- path.expand(cnaFile)
    stopifnot(file.exists(cnaFile))
  }

  if (!is.null(purityFile)){
    purityFile <- path.expand(purityFile)
    stopifnot(file.exists(purityFile))
  }

  # get vcf
  vcf <- parseVcfFile(vcfFile)

  # get cna reconstruction
  cnaRanges <- parseCnaFile(cnaFile)

  # vcaf has vcf and vaf data concatenated
  vcaf <- getVcaf(vcf, purity, cnaRanges, refGenome)
  mutTypes <- read.delim(mutTypesFile, stringsAsFactors = F)

  # strip chr if present
  mutTypes$chromosome <- unlist(strsplit(mutTypes$chromosome, "chr"))[c(F, T)]

  # match ordering
  vcaf <- vcaf[order(vcaf$pos, vcaf$chr),]
  mutTypes <- mutTypes[order(mutTypes$start, mutTypes$chromosome),]

  # sanity check
  stopifnot( all(vcaf$chr == mutTypes$chromosome) & all(vcaf$pos == mutTypes$start) )

  # set trinucs and re-order by phi
  vcaf$mutType <- mutTypes$tri
  vcaf <- vcaf[order(vcaf$phi),]

  # save intermediate is necessary
  if (saveIntermediate == TRUE){
    mut_types <- vcaf[,c("chr", "pos", "phi", "ref", "alt", "mutType")]
    write.table(mut_types, file = intermediateFile, quote = F, col.names = F, row.names = F, sep = "\t")
  }

  return( getBinCounts(vcaf, binSize, context) )


}








# [END]
