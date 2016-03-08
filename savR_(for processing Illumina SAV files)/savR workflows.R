# savR workflows
# My script to parse out UAGC desired metrics from IlluminaSAV files
# 2016 Mar 7

# savR documentation: https://github.com/bcalder/savR
# Illumina documentation: 
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/sav/sequencing-analysis-viewer-v1_8_46-guide-15066069-a.pdf
# Written using data from Illumina project "160208_SN885_0258_AHJC22BCXX"

library("savR")

# Create project
# Alternatively, setwd for wherever your project is
setwd("~/Scripts/R_scripts/savR_(for processing Illumina SAV files)/interopIlluminaData")
fc <- savR("160208_SN885_0258_AHJC22BCXX")

##################################################################################################

# Phas/Prephas(%) 
# The value used by RTA for the percentage of molecules in a cluster for which sequencing 
# falls behind (phasing) or jumps ahead (prephasing) the current cycle within a read.

# Get tile metrics
TM <- tileMetrics(fc)
# Should output a data frame with 4 columns:
  # 1) lane: Lane number
  # 2) tile: Tile ID
  # 3) code: Metrics based on the following format:
                #cluster density:	                    100
                #PF cluster density:	                101
                #number of clusters:	                102
                #number of PF clusters:	              103
                #phasing for read N:	                (200 + (N - 1) * 2), always even -> 200 (N=1) 202 (N=2) 204 (N=3)
                #prephasing for read N:	              (201 + (N - 1) * 2), always odd ->  201       203       205
                #percent aligned for read N:	        (300 + N - 1)                   ->  300       301       302
                #control lane:	                      400
  # 4) value: Value for code key

##################################################################################################

# Reads PF(M)
# The number of clusters (in millions) passing filtering.

# Sum the total pass filter number of clusters for all tiles in the lane
pfClusters(fc, 1L)

##################################################################################################

# % >= Q30
# The percentage of bases with a quality score of 30 or higher, respectively. 
# This chart is generated after the 25th cycle, and the values represent the current cycle.

# Return the ratio of clusters with a quality score less than or equal to a specified value (n) 
# for the requested lanes and cycles. In this case, n = 30.
# clusterQualityGtN(project, lane, cycle, n)
clusterQualityGtN(fc, 1L, 25L, 30L)

##################################################################################################

# Error Rate(%)
# The calculated error rate, as determined by the PhiX alignment. 
# Subsequent columns display the error rate for cycles 1–35, 1–75, and 1–100.

# Get error metrics for lane, tile and cycle. 
EM <- errorMetrics(fc)
# Should output the following:
  # lane: Lane number
  # tile: Tile ID
  # cycle: Cycle number
  # errorrate: Error rate
  # nPerfect: number of perfect reads
  # n[1-4]Error: Number of reads with 1, 2, 3 and 4 errors

