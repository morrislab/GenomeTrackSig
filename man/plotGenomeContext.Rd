% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGenomeContext.R
\name{plotGenomeContext}
\alias{plotGenomeContext}
\title{Plot the genomic profile of a tumor and other relevant genomic features.}
\usage{
plotGenomeContext(
  counts,
  trajectory,
  binSize,
  title = "",
  chr_level = F,
  cutoff = 0,
  show = T,
  chrom_states = TrackSig:::roadmap_consensus_chromstates
)
}
\arguments{
\item{counts}{un-binned dataframe of mutation counts for the sample of interest;
same dataframe as `counts` parameter in @seealso \link{GenomeTrackSig}}

\item{trajectory}{a list containing named elements "mixtures",
"changepoints", and 'binData'. See @seealso \link{GenomeTrackSig}.}

\item{binSize}{number of mutations in each bin; same value as
`binSize` parameter in @seealso \link{GenomeTrackSig}.}

\item{title}{string containing desired plot title; default is blank}

\item{chr_level}{logical whether TrackSig was run on each chromosome
separately; default is FALSE.}

\item{cutoff}{minimum proportion of bootstrap samples that must agree on a
changepoint location in order for that changepoint to be plotted; default is
0.}

\item{show}{logical whether to print the plot; default TRUE.}

\item{chrom_states}{dataframe of chromatin state information for 1Mb regions
throughout the genome.}
}
\value{
ggplot object
}
\description{
\code{plotGenomeContext} For each bin in a set of signature
mixtures, the mixture is plotted across the genome. Provided changepoints
will be highlighted. Normalized mutation density, normalized gene density, GC content,
and the distribution of consensus chromatin states (https://github.com/gerstung-lab/tensorsignatures) are also plotted at each bin.
}
