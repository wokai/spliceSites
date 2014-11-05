
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  File   : allClasses.r                                                        #
#  Date   : 05.Nov.2012                                                         #
#  Content: Class definitions for package spliceSites                           #
#  Version: 0.7.0                                                               #
#  Author : W. Kaisers                                                          #
#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Changelog
#  03.Jul.13  Added unloading of spliceSites.so
#  09.Jul.13  Added maxEnt class
#  10.Jul.13  Removed gapRanges, gapProbes, and derived (aaX,dnaX) classes.
#  24.Jul.13  Changed alingGaps to gapSites
# 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

.onUnload<-function(libpath) { library.dynam.unload("spliceSites",libpath) }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  Declarations for cRanges                                                                         #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Centered Range: contains id and position
setClass("cRanges",representation(dt="data.frame"))
# Centered Range with DNAStringSet
setClass("cdRanges",representation(seq="DNAStringSet"),contains="cRanges")
# Centered Range with AAStringSet
setClass("caRanges",representation(seq="AAStringSet"),contains="cRanges")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  Declarations for x-Ranges                                                                        #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Model:                                                          +
#                                                                 +
#                                                                 +
#        1000 > > >  ascending reference positions > > > 2000     +
#                                                                 +
#                 left                        right               +
#                                                                 +
#             lstart     lend            rstart     rend          +
#             ###############<- gaplen ->###############          +
#             |    lfeat    |    gap     |    rfeat    |          +
#                                                                 +
#                                                                 +
# alt_left :             alt             group                    +
# alt_right:           group             alt                      +
#                                                                 +
#                                                                 +
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


setClassUnion("dataFrameOrNULL",c("data.frame","NULL"))

# Additionally contains information about Align-numbers plus eventually annotation data.
setClass("gapSites",representation(dt="data.frame",nAligns="numeric",nAlignGaps="numeric",lann="dataFrameOrNULL",rann="dataFrameOrNULL",profile="dataFrameOrNULL")
         ,prototype(dt=data.frame(),nAligns=0,nAlignGaps=0,lann=NULL,rann=NULL,profile=NULL))
# Additionally stores probe and Annotation overlap information in dt
setClass("annGapSites",contains="gapSites")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  Declarations for dnaGapSites                                                                     #
#  Additionally contain DNAStringSet                                                                #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("dnaGapSites",representation(seq="DNAStringSet"),contains="gapSites")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  Declarations for aaGapSites                                                                      #
#  Additionally contain AAStringSet                                                                 #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("aaGapSites",representation(seq="AAStringSet"),contains="gapSites")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  keyProfiler:                                                                                     #
#  counts occurrence of profile (prof) items defined by factor levels in profile table              #
#  in successively added tables by keys defined in object@ev$dt                                     #
#  unique: allows adding keys only once for each profile line                                       #
#  values: when given, conditions are not simply counted but the given values                       #
#          are summed up for each condition.                                                        #
#          values must either be absent or present for all accumulated data                         #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("keyProfiler",representation(ev="environment",unique="logical",counted="logical",useValues="logical"))


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  SpliceCountSet                                                                                   #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("SpliceCountSet",contains="ExpressionSet")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  maxEnt                                                                                           #
#  Calculates splice site scores                                                                    #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("maxEnt",
  representation(
    ev="environment",
    basedir="character"),
  prototype=prototype(
     ev=new.env(),
     basedir="."
))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  hbond                                                                                            #
#  Calculates H-bond scores                                                                         #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("hbond",
         representation(
           ev="environment",
           basedir="character"),
         prototype=prototype(
           ev=new.env(),
           basedir="."
         ))


