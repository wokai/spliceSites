
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  c for gapRanges
setMethod("c","cRanges",function(x, ..., recursive=FALSE)
{
  if (recursive)
    stop("'recursive' mode not supported!")
  cr<-new(class(x))
  args <- unname(list(x, ...))
  cr@dt<-do.call(rbind,lapply(args,getDataFrame))
  if(.hasSlot(x,"seq"))
    cr@seq<-do.call(c,lapply(args,getSequence))  
  return(cr)
})

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  c for gapSites

setMethod("c","gapSites",function(x,...,recursive=FALSE)
{
  if(recursive)
    warning("recursive mode not supported.")
  
  ga<-new(class(x))
  l<-list(x,...)
  #nr<-sum(unlist(lapply(l,nrow)))
  args<-unname(l)
  ga@dt<-do.call(rbind,lapply(args,getDataFrame))
  
  
  ga@lann<-do.call(rbind,lapply(args,getLann))
  ga@rann<-do.call(rbind,lapply(args,getRann))
  ga@nAligns<-sum(unlist(lapply(args,nAligns)))
  ga@nAlignGaps<-sum(unlist(lapply(args,nAlignGaps)))
  
  # Recalculate gptm and rpmg
  gptm_fac<-1e7/ga@nAligns
  rpmg_fac<-1e6/ga@nAlignGaps
  ga@dt$id<-1:nrow(ga@dt)
  ga@dt$gptm<-round(ga@nAligns*gptm_fac,3)
  ga@dt$rpmg<-round(ga@nAligns*rpmg_fac,3)  
  
  # Derived: dnaGapSites, aaGapSites
  if(.hasSlot(x,"seq"))
    ga@seq<-do.call(c,lapply(args,getSequence))  
  
  return(ga)
})