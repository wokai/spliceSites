
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## c-methods
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##  c for gapRanges
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

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

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##  c for gapSites
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("c","gapSites",function(x,...,recursive=FALSE)
{
    if(recursive)
        warning("recursive mode not supported.")
    
    ga<-new(class(x))
    
    l<-list(x,...)
    args<-unname(l)
    
    # rbind position and annotation tables
    getAnn <- function(object) {return(object@annotation)}
    getDataFrame <- function(object) {return(object@dt)}
    
    ga@dt<-do.call(rbind,lapply(args,getDataFrame))
    ga@annotation<-do.call(rbind,lapply(args,getAnn))
    
    # sum up nAligns and nAlignGaps
    ga@nAligns<-sum(unlist(lapply(args,nAligns)))
    ga@nAlignGaps<-sum(unlist(lapply(args,nAlignGaps)))
    
    
    ## Recalculate gptm and rpmg
    gptm_fac<-1e7/ga@nAligns
    rpmg_fac<-1e6/ga@nAlignGaps
    ga@dt$id<-1:nrow(ga@dt)
    ga@dt$gptm<-round(ga@nAligns*gptm_fac,3)
    ga@dt$rpmg<-round(ga@nAligns*rpmg_fac,3)  
    
    
    ## Derived: dnaGapSites, aaGapSites
    if(.hasSlot(x,"seq"))
        ga@seq<-do.call(c,lapply(args,getSequence))  
    
    return(ga)
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
