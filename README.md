# GenomeTrackSig R package (dev)
Morris Lab, Sloan Kettering Institute. R package for GenomeTrackSig. 

To cite, please see https://www.biorxiv.org/content/10.1101/2022.01.23.477261v1

# Vignette 
Coming soon

# Dependencies 
R >= 3.3.3

This package imports the following other R packages:

+ reshape2 >= 1.4.3 (CRAN)
+ ggplot2 >= 3.2.0 (CRAN)
+ NMF >= 0.21.0 (CRAN)
+ assertthat >= 0.2.1 (CRAN)
+ BSgenome.Hsapiens.UCSC.hg19 >= 1.4.0 (Bioconductor)
+ GenomicRanges >= 1.26.4 (Bioconductor)
+ Biostrings >= 2.42.1 (Bioconductor)
+ SummarizedExperiment >= 1.4.0 (Bioconductor)
+ VariantAnnotation >= 1.20.3 (Bioconductor)
+ grid >= 3.3.3 (CRAN)
+ progress >= 1.2.2 (CRAN)

# Load the package in R
`devtools::install_github("morrislab/GenomeTrackSig")`

# Demo
Using the example data provided in `extdata/`, the following code will plot the signature profile, and return the fitted mixture of signatures for each bin, the bins where changepoints were detected, and the ggplot object.

1. First, restrict the list of signatures to fit exposure for. This is recommended for improving speed by making the model smaller. Here, we choose a threshold of 5%, meaning that signatures with exposure under this across all timepoints will not be fit. 

```r
library(GenomeTrackSig)
library(ggplot2)

counts = system.file(package = "GenomeTrackSig", "extdata/Example_counts.csv")

detectedSigs <- detectActiveSignatures(counts, binSize = 200)
```
2. Next, we compute the signature activity profile across all regions of the genome. 
```r
set.seed(1224)

traj <- GenomeTrackSig(counts = counts, sampleID = "example", activeInSample =     
                      detectedSigs, binSize = 200)
```

3. Plot the activity profile. 

```r
plotGenomeProfile(traj, chr_level=F, cutoff=0)
```


# Note

Some users may have plotting issues with `GenomeTrackSig` if `ggplot2` is not explicitly loaded with `library(ggplot2)`. We are experiencing a bug that has been [previously described](https://github.com/tidyverse/ggplot2/issues/663) for `ggplot2`.
