# savR practice
# Used for Illumina SAV files from project "160208_SN885_0258_AHJC22BCXX"
# 2016-Feb-18

# Reference materials:
# https://www.bioconductor.org/packages/release/bioc/manuals/savR/man/savR.pdf
# https://www.bioconductor.org/packages/release/bioc/vignettes/savR/inst/doc/savR.pdf

library("savR")

# Create project
setwd("~/Documents/interopIlluminaData") #Alternatively, setwd for wherever your project is
fc <- savR("160208_SN885_0258_AHJC22BCXX")

# Return the ratio of clusters with a quality score less than or equal to a specified value (n) 
# for the requested lanes and cycles
# clusterQualityGtN(project, lane, cycle, n)

clusterQualityGtN(fc, 1L, 25L, 30L)

# Sum the total number of clusters for all tiles in the lane
clusters(fc, 1L)

# Sum the total pass filter number of clusters for all tiles in the lane
pfClusters(fc, 1L)

# Generate a boxplot of the numbers of clusters and the number of Illumina pass-filter clusters 
# per tile and lane
pfBoxplot(fc)

# Get number of cycles
cycles(fc)

# Get the number of sequencing reads (excluding index reads)
directions(fc)

# Get error metrics for lane, tile and cycle. 
EM <- errorMetrics(fc)
# Should output the following:
  # lane: Lane number
  # tile: Tile ID
  # cycle: Cycle number
  # errorrate: Error rate
  # nPerfect: number of perfect reads
  # n[1-4]Error: Number of reads with 1, 2, 3 and 4 errors

# Accessor to obtain information about the characteristics of the flowcell from an Illumina sequencing run
flowcellLayout(fc)
# Should output the following:
  # lanecount: Number of lanes on the flowcell
  # surfacecount: Number of surfaces
  # swathcount: Number of imaging swaths
  # tilecount: Number of tiles per swath
  # sectionperlane: Number of sections per lane (NextSeq)
  # lanepersection: Number of lanes per section (NextSeq) tilenamingconvention: Description of deviation from original formatting layout

# Accessor to obtain information about the reads of a particular Illumina sequencing run
reads(fc)
# Should output:
  # number: the index of this read in sequencing
  # cycles: number of cycles in this read
  # index: logical representing whether or not this read is an index read

# Generate a plot for a given cycle of the percentage of clusters in each tile that are >= Q30 (Currently not working)
# Getting the following error message:
# "trying to get slot "data" from an object of a basic class ("NULL") with no slots"
plotQGT30(fc, cycle = 1L)

# Get tile metrics
TM <- tileMetrics(fc)
# Should output:
  # lane: Lane number
  # tile: Tile ID
  # code: Metrics based on the following format:
      #cluster density:	                    100
      #PF cluster density:	                101
      #number of clusters:	                102
      #number of PF clusters:	              103
      #phasing for read N:	                (200 + (N - 1) * 2), always even -> 200 (N=1) 202 (N=2) 204 (N=3)
      #prephasing for read N:	              (201 + (N - 1) * 2), always odd ->  201       203       205
      #percent aligned for read N:	        (300 + N - 1)                   ->  300       301       302
      #control lane:	                      400
  # value: Value for code key

