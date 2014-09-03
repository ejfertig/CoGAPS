# GAPS: retired function, originally for link to JAGS

#'\code{GAPS} is retired as of version 2, use gapsRun instead
#'@export

GAPS <- function(data, unc,
                 outputDir, outputBase="",
                 sep = "\t", isPercentError=FALSE,
                 numPatterns,
                 MaxAtomsA=2^32, alphaA=0.01,
                 MaxAtomsP=2^32, alphaP=0.01,
                 SAIter=1000000000, iter = 500000000, thin=-1,
                 verbose=TRUE, keepChain = FALSE) {
  
    print("Please Use gapsRun Function in This Version")
    print("of CoGAPS.  See Help Files for New Call")
}
                          
