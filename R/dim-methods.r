
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## dim-methods
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("dim","cRanges",function(x)  dim(x@dt))
setMethod("dim","gapSites",function(x) dim(x@dt))


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
