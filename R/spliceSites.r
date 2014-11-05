
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
##                                                                                              ##
##  Project   :   spliceSites                                                                   ##
##  Created   :   29.Okt.2012                                                                   ##
##  Author    :   W. Kaisers                                                                    ##
##  Content   :   Funktionality for working with RNA-seq splice-sites data                      ##
##                based von rbamtools and refGenome packages.                                   ##
##                for usage in R                                                                ##
##  Version   :   1.1.1                                                                         ##
##                                                                                              ##
##  Changelog :                                                                                 ##
##  12.11.12  :   Corrected false calculation of gaplen in merge.alignGaps                      ##
##  12.12.12  :   Added Dependency on BiocGenerics                                              ##
##  03.12.12  :   Corrected and checked Proteomic pipeline: lJunc,rJunc,lCodons,rCodons,        ##
##                lrJunc,lrCodons                                                               ##
##  04.12.12  :   Added keyProfiler class                                                       ##
##  11.01.13  :   0.1.1                                                                         ##
##  06.02.13  :   Added readMergedBamSites function which uses rbamtools:::bamSiteLists         ##
##  18.02.13  :   0.2.1                                                                         ##
##  28.05.13  :   0.3.1                                                                         ##
##  03.06.13  :   0.3.3 added suppression of scientific notation in write.annDNAtables          ##
##  03.06.13  :   0.3.4 readExpSet                                                              ##
##  03.07.13  :   0.4.0 Added truncateSeq (C-routine) and changed trypsinCleave                 ##
##                to C-implementation                                                           ##
##  09.07.13  :   0.5.0 Added maxEnts                                                           ##
##  10.07.13  :   0.5.1 Changed position entries from 0-based to 1-based                        ##
##  10.07.13  :   0.5.2 Removed gapRanges, gapProbes and derived (aaX, dnaX) classes            ##
##  18.07.13  :   0.6.0 Renamed gapAligns to alignGaps, removed readMergedBamGapProbes          ##
##  24.07.13  :   0.7.0 Renamed alignGaps to gapSites                                           ##
##  31.07.13  :   0.8.0 Removed dependency on data.table, C-version for alt_group               ##
##  09.08.13  :   0.99.0 All C routines valgrind checked                                        ##
##  21.08.13  :   0.99.7 Review by Marc Carlson:                                                ##
##                          Added biocViews in DESCRIPTION                                      ##
##                          (RNAseq,GeneExpression,DifferentialExpression,Proteomics)           ##
##                          Added Collate entry in DESCRIPTION                                  ##
##                          Added R-registration of C-routines                                  ##
##                          Switched calloc to R_alloc in C (valgrind tested)                   ##
##  21.08.13  :   0.99.8  Added functions for hbond score (valgrind tested)                     ##
##  23.08.13  :   0.99.9  Submission update for Bioconductor                                    ##
##  23.08.13  :   0.99.10 Added get_dna_nmers (valgrind tested)                                 ##
##  03.12.13  :   1.1.1  Changed annotate.gapSites function                                     ##
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##

.onUnload<-function(libpath) { library.dynam.unload("spliceSites",libpath) }


strandlevels<-c("+","-","*")

# Unexported function for creating index values on alphabet (->readExpSet)
alphabetIndex<-function(idxLen,alpha)
{return(.Call("get_alph_index",as.integer(idxLen),alpha,PACKAGE="spliceSites"))}

# 36 character alphabet
index_alphabet<-paste(paste(0:9,collapse=""),paste(letters,collapse=""),sep="")


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#                                                                                                   #
#  Declarations for cRanges:                                                                        #
#  Centered Ranges                                                                                  #
#  Ranges which contain a pointer to a "position" of interest inside the range.                     #
#  Intended to represent a range around an exon-intron boundary                                     #
#  where "position" represents the                                                                  #
#  1-based position of the last exon nucleotide                                                     #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("initialize","cRanges",function(.Object,
                                          seqid=NULL,
                                          start=NULL,end=NULL,width=NULL,
                                          strand=NULL,
                                          position=0L,id=NULL)
{
  # Allows to create an empty object
  if(is.null(start) && is.null(end))
    return(.Object)
  
  n<-length(start)
  if(is.null(id))
    id<-1:n
  if(!is.null(width))
    end<-start+width-1L
  
  # +++++++++++++++++++++++++++++++++++
  #  strand
  if(is.null(strand))
    strand<-factor(rep("*",n),levels=strandlevels)
  if(length(strand)==1)
    strand=factor(rep(strand,n),levels=strandlevels)
  
  # +++++++++++++++++++++++++++++++++++
  if(length(seqid)==1)
    seqid=factor(rep(seqid,n))
  seqid=factor(seqid)
  
  dfr<-data.frame(seqid=seqid,start=start,end=end,strand=strand,position=as.integer(position),id=id)
  .Object@dt<-dfr[order(dfr$seqid,dfr$start,dfr$end),]
  return(.Object)
})

setMethod("initialize","cdRanges",function(.Object,
                                           seqid=NULL,
                                           start=NULL,end=NULL,width=NULL,
                                           strand=NULL,
                                           position=0L,id=NULL,seq=NULL)
{
  .Object<-callNextMethod(.Object,seqid,start,end,width,strand,position,id)
  n<-nrow(.Object@dt)
  if(is.null(seq))
    .Object@seq<-DNAStringSet()
  else if(is.character(seq))
    .Object@seq<-DNAStringSet(seq)
  else if(is(seq,"DNAStringSet"))
    .Object@seq<-seq
  else
    stop("[initialize.cdRanges] seq must be character or DNAStringSet!") 
  return(.Object)
})

setValidity("cRanges",
            function(object){
              if(any(object@dt$end-object@dt$start<0))
                return("end must be >= start!")
              
              # + + + + + + + + + + + + 
              # Range:     1 2 3 4
              # Position:      3
              # End-Start=3
              # + + + + + + + + + + + +
              
              if(any(object@dt$position>(object@dt$end-object@dt$start+1)))
                return("position must be inside of range (i.e. <=end-start+1)!")
              else
                return(TRUE)
            })


setMethod("seqid", "cRanges",function(x)return(x@dt$seqid))
setMethod("start", "cRanges",function(x) return(x@dt$start))
setMethod("end",   "cRanges",function(x)return(x@dt$end))
setMethod("id",    "cRanges",function(x)return(x@dt$id))
setMethod("strand","cRanges",function(x,...)return(x@dt$strand))
setMethod("count", "cRanges",function(x)return(nrow(x@dt)))
setMethod("width", "cRanges",function(x)return(x@dt$end-x@dt$start+1L))
setMethod("getSequence","cdRanges",function(x) return(x@seq))
setMethod("getSequence","caRanges",function(x) return(x@seq))

setMethod("sortTable","cRanges",function(x){
  o<-order(x@dt$seqid,x@dt$start,x@dt$end)
  cr<-new(class(x))
  cr@dt<-x@dt[o,]
  if(.hasSlot(x,"seq"))
    cr@seq<-x@seq[o]
  return(cr)
})           

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# xCodons functions
# Do upstream frame-shift and downstream trim to full codon.
# Position values are updated
# (position = 1-based position of last exon nucleotide)
# Eventually overwrite strand column
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
# lCodons: 
# Left side frame-shift and right-trim to full codon
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("lCodons","cRanges",function(x,frame=1,keepStrand=TRUE){
  if(!is.element(frame,1:3))
    stop("[lCodons.cRanges] frame must be one of 1,2,3!")
  
  # Convert frame into 0-based shift value
  nframe<-as.integer(frame-1)

  # Create return object
  cr<-new("cRanges")
  cr@dt<-x@dt

  cr@dt$start<-cr@dt$start+nframe
  
  # Calculate width and truncate
  # width to full codon length
  w<-cr@dt$end-cr@dt$start+1L
  w<-w-w%%3
  
  cr@dt$end<-cr@dt$start+w-1L
  cr@dt$position<-cr@dt$position-nframe
  
  n<-nrow(cr@dt)
  cr@dt$frame<-rep(frame,n)
  if(!keepStrand)
    cr@dt$strand<-factor(rep("+",n),levels=strandlevels)
  
  if(any(cr@dt$position<0))
    message("[lCodons.cRanges] Found position <0!")
  if(any(cr@dt$position>w))
    message("[lCodons.cRanges] Found position > width!")
  return(cr)
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# 
# rCodons: 
# Right side frame-shift and left-trim to full codon
# Intended to be used in combination with reverseComplement (later on)
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("rCodons","cRanges",function(x,frame=1,keepStrand=TRUE){
  if(!is.element(frame,1:3))
    stop("[rCodons.cRanges] frame must be one of 1,2,3!")
  
  # Convert frame into 0-based shift value
  nframe<-as.integer(frame-1)

  # Create return object
  cr<-new("cRanges")
  cr@dt<-x@dt

  # Right shift
  cr@dt$end<-cr@dt$end-nframe
  
  # Calculate width and truncate
  # width to full codon length
  w<-cr@dt$end-cr@dt$start+1L
  w<-w-w%%3
    
  # Left truncation
  cr@dt$start<-cr@dt$end-w+1L
  # Correct position
  cr@dt$position<-cr@dt$position-nframe
  
  n<-nrow(cr@dt)  
  cr@dt$frame<-rep(frame,n)
  if(!keepStrand)
    cr@dt$strand<-factor(rep("-",n),levels=strandlevels) 
  
  if(any(cr@dt$position<0))
    message("[rCodons.cRanges] Found position <0!")
  if(any(cr@dt$position>w))
    message("[rCodons.cRanges] Found position > width!")
  return(cr)
})


setMethod("dnaRanges","cRanges",function(x,dnaset,useStrand=TRUE,removeUnknownStrand=TRUE,verbose=TRUE...){
  if(!is(dnaset,"DNAStringSet"))
    stop("[dnaRanges.cRanges] dnaset must be DNAStringSet!")
  if(!is.logical(useStrand))
    stop("[dnaRanges.cRanges] useStrand must be logical!")
  
  dnames<-names(dnaset)
  if(is.null(dnames))
    stop("[dnaRanges.cRanges] names(dnaset) must not be NULL!")

  bm<-Sys.localeconv()[7]
  if(removeUnknownStrand&&useStrand)
  {
    lstr<-x@dt$strand!="*"
    nlstr<-sum(lstr)
    if(nlstr==0)
      stop("[dnaRanges.cRanges] Removing all(!) rows from cRanges due to unknown strand!")
    if(verbose)
      message("[dnaRanges.cRanges] Removing ",format(length(lstr)-nlstr,big.mark=bm),
              " rows from cRanges due to unknown strand (rest: ",format(nlstr,big.mark=bm)," rows).")
    x@dt<-x@dt[lstr,]
  }
  
  n<-length(dnames)
  # New table as tamplate
  ans<-x@dt[0,]
  ans$seq<-character()
  for(i in 1:n)
  {
    lg<-as.logical(seqid(x)==dnames[i])
    if(any(lg))
    {
      dt<-x@dt[x@dt$seqid==dnames[i],]
      dt$seq<-as.character(Views(dnaset[[i]],dt$start,dt$end))
      # Add data for seqid to data.frame
      ans<-rbind(dt,ans)
    }
  }
  
  # Sort ans
  ans<-ans[order(ans$seqid,ans$start,ans$end),]
  if(nrow(ans)==0)
    stop("[dnaRanges.cRanges] No match of sequence names between seqid(cRanges) and names(dnaset)!")

  dss<-DNAStringSet(ans$seq)
  ans$seq<-NULL
  if(useStrand)
  {
    lg<-as.logical(ans$strand=="-")
    dss[lg]<-reverseComplement(dss[lg])
    if(verbose)
    {
      if(removeUnknownStrand)
        message("[dnaRanges.cRanges] Reversed ",format(sum(lg),big.mark=bm)," sequences.")      
      else
        message("[dnaRanges.cRanges] useStrand: '*' is treated as '+'. Reversed ",
                format(sum(lg),big.mark=bm)," sequences.")
    }
  }
  res=new("cdRanges")
  res@dt<-ans
  res@seq<-dss
  return(res)
})


setMethod("extractRange","gapSites",function(object,seqid,start,end){
  if(!is.character(seqid))
    stop("[extractRange.gapSites] seqid must be character!")
  if(!is.numeric(start))
    stop("[extractRange.gapSites] start must be numeric!")
  start<-as.integer(start)
  if(!is.numeric(end))
    stop("[extractRange.gapSites] end must be numeric!")
  end<-as.integer(end)
  if(start>=end)
    stop("[extractRange.gapSites] end must be greater than start!")
  
  sq<-object@dt$seqid==seqid
  st<-object@dt$lstart>=start
  se<-object@dt$rend<=end
  res<-sq&st&se
  
  exr<-new(class(object))
  exr@dt<-object@dt[res,]
  exr@nAligns<-object@nAligns
  exr@nAlignGaps<-object@nAlignGaps
  
  mtc<-match(exr@dt$id,object@rann$id)
  exr@rann<-object@rann[mtc,]
  mtc<-match(exr@dt$id,object@lann$id)
  exr@lann<-object@lann[mtc,]
  
  if(.hasSlot(object,"seq"))
    exr@seq<-object@seq[res]
  return(exr)
})

setMethod("extractRange","cRanges",function(object,seqid,start,end)
{
  if(!is.character(seqid))
    stop("[extractRange.cRanges] seqid must be character!")
  if(!is.numeric(start))
    stop("[extractRange.cRanges] start must be numeric!")
  start<-as.integer(start)
  if(!is.numeric(end))
    stop("[extractRange.cRanges] end must be numeric!")
  end<-as.integer(end)
  if(start>=end)
    stop("[extractRange.cRanges] end must be greater than start!")
  
  sq<-object@dt$seqid==seqid
  st<-object@dt$start>=start
  se<-object@dt$end<=end
  res<-sq&st&se
  exr<-new(class(object))
  exr@dt<-object@dt[res,]
  # cdRanges and caRanges:
  if(.hasSlot(object,"seq"))
    exr@seq<-object@seq[res]
  return(exr)
})


setMethod("extractByGeneName","cRanges",function(object,geneNames,src,...)
{
  if(missing(geneNames))
    stop("[extractByGeneName.cRanges] geneNames argument is not optional!")
  if(missing(src))
    stop("[extractByGeneName.cRanges] src argument is not optional!")
  if(!is.character(geneNames))
    stop("[extractByGeneName.cRanges] geneNames must be character!")  
  if(!is(src,"refGenome"))
    stop("[extractByGeneName.cRanges] src must be of class 'refGenome'!")
  
  reg<-extractByGeneName(src,geneNames)
  gpos<-getGenePositions(reg)
  n<-nrow(gpos)
  if(n==0)
    stop("[extractByGeneName.cRanges] Empty gene position table (no matches in geneNames?)!")
  res<-extractRange(object,seqid=as.character(gpos$seqid[1]),start=gpos$start[1],end=gpos$end[1])
  if(n>1)
  {
    for(i in 1:n)
      res<-c(res,extractRange(object,seqid=as.character(gpos$seqid[i]),start=gpos$start[i],end=gpos$end[i]))
  }
  return(res)
})

setMethod("extractByGeneName","gapSites",function(object,geneNames,src,...){
  if(missing(geneNames))
    stop("[extractByGeneName.gapSites] geneNames argument is not optional!")
  if(missing(src))
    stop("[extractByGeneName.gapSites] src argument is not optional!")
  if(!is.character(geneNames))
    stop("[extractByGeneName.gapSites] geneNames must be character!")  
  if(!is(src,"refGenome"))
    stop("[extractByGeneName.gapSites] src must be of class 'refGenome'!") 
  
  reg<-extractByGeneName(src,geneNames)
  gpos<-getGenePositions(reg)
  n<-nrow(gpos)
  if(n==0)
    stop("[extractByGeneName.gapSites] Empty gene position table (no matches in geneNames?)!")
  
  res<-extractRange(object,seqid=as.character(gpos$seqid[1]),start=gpos$start[1],end=gpos$end[1])
  if(n>1)
  {
    for(i in 1:n)
      res<-c(res,extractRange(object,seqid=as.character(gpos$seqid[i]),start=gpos$start[i],end=gpos$end[i]))
  }
  return(res)
})


setMethod("seqlogo","cdRanges",function(x,strand="+",useStrand=TRUE,...){
  if(length(x@seq)==0)
    stop("[seqlogo.cdRanges] length(DNAStringSet)=0!")
  
  # Generates seqLogo for whole DNAStringSet
  bm<-Sys.localeconv()[7]
  if(useStrand)
  {
    strandseqs<-x@dt$strand==strand
    message("[seqlogo.cdRanges] +: ",format(sum(strandseqs),big.mark=bm),", total: ",format(nrow(x@dt),big.mark=bm),".")
    if(sum(strandseqs)==0)
      stop("[seqlogo.cdRanges] No range found for strand '",strand,"'!")
    cs<-consensusMatrix(x@seq[strandseqs],as.prob=T,baseOnly=TRUE)
  }
  else
    cs<-consensusMatrix(x@seq,as.prob=T,baseOnly=TRUE)
  pwm<-makePWM(cs[1:4,])
  seqLogo(pwm,ic.scale=F)    
})

setMethod("translate","cdRanges",function(x){
  
  # Exclude data with N's in Sequence
  af<-alphabetFrequency(x@seq)
  withN<-af[,"N"]>0
  noN<-!withN
  bm<-Sys.localeconv()[7]
  if(sum(withN)>0)
  {
    message("[translate.cdRanges] Excluding ",format(sum(withN),big.mark=bm),"/",
        format(nrow(x@dt),big.mark=bm)," rows because of N's!")
    x@dt<-x@dt[noN,]
    x@seq<-x@seq[noN]
  }
  
  ca<-new("caRanges")
  ca@dt<-x@dt
  
  # + + + + + + + + + + + + + + + + + + #
  # Range:    1 2 3 | 4 5 6 | 7 8 9
  # Position:         4
  # Corrected position: 2 (4L %/% 3L = 1L)
  # + + + + + + + + + + + + + + + + + + #
  
  ca@dt$position<-(ca@dt$position %/% 3L)+1L
  ca@seq<-translate(x@seq)
  return(ca)
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# truncSeq - functions
# Truncation functions for amino-acid containing objects
# Intended to cut sequences at stop-codons and recalculate dependent positions
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

truncate_seq<-function(seq,pos,id,rme=TRUE,trunc=42L)
{
  if(!is.character(seq))
    stop("[truncate_seq] seq must be character!")
  if(!is.numeric(pos))
    stop("[truncate_seq] pos must be numeric!")
  if(length(seq)!=length(pos))
    stop("[truncate_seq] seq and pos must have same length!")

  if(missing(id))
    id<-1:length(seq)
  else {
    if(!is.numeric(id))
      stop("[truncate_seq] id must be numeric!")
    id<-as.integer(id)
  }

  pos<-as.integer(pos)
  if(!is.logical(rme))
    stop("[truncate_seq] rme must be logical!")
  if(!is.numeric(trunc))
    stop("[truncate_seq] trunc must be numeric!")
  trunc<-as.integer(trunc)
  
  ret<-.Call("trunc_pos",id,pos,seq,trunc[1],PACKAGE="spliceSites")
  if(rme)
    ret<-ret[ret$lseq>0,]
  return(ret)  
}

trnctsq<-function(x,rme=TRUE,trunc=42L){
  # rme: remove empty seqs
  if(!is.logical(rme))
    stop("[truncateSeq] rme must be logical!")
  
  ret<-.Call("trunc_pos",x@dt$id,as.integer(x@dt$position),as.character(x@seq),as.integer(trunc[1]),PACKAGE="spliceSites") 
  rr<-new(class(x))
  
  if(rme[1])
  {
    nn<-ret$lseq>0
    rr@dt<-x@dt[nn,]
    rr@dt$position<-ret$pos[nn]
    rr@dt$lseq<-ret$lseq[nn]
    rr@seq<-AAStringSet(ret$seq[nn])
    return(rr)
  }
  
  rr@dt<-x@dt
  rr@dt$position<-ret$pos
  rr@dt$lseq<-ret$lseq
  rr@seq<-AAStringSet(ret$seq)
  return(rr)  
}

setMethod("truncateSeq","caRanges",   function(x,rme=TRUE,trunc=42L){trnctsq(x,rme=TRUE,trunc=trunc)})
setMethod("truncateSeq","aaGapSites",function(x,rme=TRUE,trunc=42L){trnctsq(x,rme=TRUE,trunc=trunc)})
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# trypsinCleave: Performs in silico trypsinization
# The C-routine implements the "Keil"-rule, where sites are described by
# the regex  "[RK](?!P)". The cut position is between [RK] and
# the following character.
# 
# The function returns the fragment which contains the position depicted
# exon-intron boundary.
#
# Sequences where the trypsinization product is shorter then minLen are
# excluded from the result.
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

trclv<-function(x,minLen=5,...){
  # position: 1-based position of last AA on left side
  # position must be given as 0-based for "silic_tryp_pos"
  
  ret<-.Call("silic_tryp_pos",as.integer(x@dt$id),as.integer(x@dt$position-1),as.character(x@seq),PACKAGE="spliceSites")
  nn<-ret$lseq>=minLen  
  
  rr<-new(class(x))
  rr@dt<-x@dt[nn,]
  rr@dt$position<-ret$pos[nn]
  rr@dt$lseq<-ret$lseq[nn]
  rr@seq<-AAStringSet(ret$seq[nn])
  
  # Remove too short AA sequences
  bm<-Sys.localeconv()[7]
  if(sum(nn)==0)
    warning("[trypsinCleave] No AA-Sequence left by filter width >= minLen (=",minLen,")!")
  message("[trypsinCleave] Removing ",format(sum(!nn),big.mark=bm)," sequences shorter than ",minLen,".")
  return(rr)  
}

setMethod("trypsinCleave","caRanges",   function(x,minLen=5,...){trclv(x,minLen,...)})
setMethod("trypsinCleave","aaGapSites",function(x,minLen=5,...){trclv(x,minLen,...)})


silic_tryp<-function(seq,pos,id)
{
  if(!is.character(seq))
    stop("[silic_tryp] seq must be character!")
  if(!is.numeric(pos))
    stop("[silic_tryp] pos must be numeric!")
  if(length(seq)!=length(pos))
    stop("[silic_tryp] seq and pos must have same length!")

  if(missing(id))
    id<-1:length(seq)
  else {
    if(!is.numeric(id))
      stop("[silic_tryp] id must be numeric!")
    id<-as.integer(id)   
  }
  pos<-as.integer(pos)
  return(.Call("silic_tryp_pos",id,pos,seq,PACKAGE="spliceSites"))
}

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("write.files","caRanges",function(x,path,filename,...){
  if(!file.exists(path))
    dir.create(path)
  tablefile<-file.path(path,paste(filename,".csv",sep=""))
  fastafile<-file.path(path,paste(filename,".fa",sep=""))
  x@dt$fasta_id<-1:nrow(x@dt)
  x@dt$seq<-as.character(x@seq)
  names(x@seq)<-paste(x@dt$fasta_id,filename,sep="|")
  write.table(x@dt,tablefile,sep=";",dec=",",row.names=F,...)
  writeXStringSet(x@seq,fastafile)
  return(invisible())
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#                                                                                                   #
#  Definitions for gapSites class                                                                  #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

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

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
#  Constructor function:
#  Makes 'gapSites' object out of raw coordinate values.
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

gapSites<-function(seqid=factor(),lstart=integer(),lend=integer(),
                    rstart=integer(),rend=integer(),gaplen,strand,
                    nr_aligns=1,nAligns=sum(nr_aligns),
                    nAlignGaps=sum(nr_aligns),nProbes=1)
{
  n<-length(seqid)
  if(n==0)
    return(new("gapSites"))
  if(!is(seqid,"factor"))
    seqid<-factor(seqid)  
  
  if(length(lstart)!=n)
    stop("[gapSites] length(lStart) ", length(lstart), " unequal length(seqid) ",n,"!")
  if(length(lend)!=n)
    stop("[gapSites] length(lend) ",   length(lend),   " unequal length(seqid) ",n,"!")
  if(length(rstart)!=n)
    stop("[gapSites] length(rstart) ", length(rstart), " unequal length(seqid) ",n,"!")
  if(length(lend)!=n)
    stop("[gapSites] length(lend) ",   length(lend),   " unequal length(seqid) ",n,"!")  
  if(missing(gaplen))
    gaplen<-rstart-lend-1L
  if(length(gaplen)!=n)
    stop("[gapSites] length(gaplen) ",length(gaplen)," unequal length(seqid) ",n,"!")
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  #  Align number values
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  if(length(nr_aligns)==1)
    nr_aligns<-rep(nr_aligns,n)
  if(length(nr_aligns)!=n)
    stop("[gapSites] length(nr_aligns) ",length(nr_aligns), " unequal length(seqid) ",n,"!")
  
  if(length(nAligns)!=1)
    stop("[gapSites] length(nAligns) ",length(nAligns)," must be 1!")
  if(length(nAlignGaps)!=1)
    stop("[gapSites] length(nAlignGaps) ",length(nAligns)," must be 1!")
  
  nAligns<-as.double(nAligns)
  nAlignGaps<-as.double(nAlignGaps)
  nProbes<-as.double(nProbes)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  #  strand
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  if(missing(strand))
    strand<-factor(rep("*",n),levels=strandlevels)  
  if(length(strand)==1)
  {
    if(strand==1 || strand=="-")
      strand<-factor(rep("-",n),levels=strandlevels)
    else
      strand<-factor(rep("+",n),levels=strandlevels)
  }
  else if(length(strand)==n)
  {
    lgl<-logical(n)
    if(is.character(strand))
      strand<-factor(ifelse(strand=="+","+","-"),levels=strandlevels)
    else if(is.numeric(strand))
      strand<-factor(ifelse(strand==0,"+","-"),levels=strandlevels)
    else if(!is(strand,"factor"))
    {
      if(!setequal(levels(strand),strandlevels))
        stop("[gapSites] strand levels must be '+-*'!")
    }
  }
  else
    stop("[gapSites] length(strand) must be 1 or length(seqid)=",length(seqid),"!")
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  #  rpmg and gptm (align number derived values)
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  rpmg_fac<-1e6/nAlignGaps
  gptm_fac<-1e7/nAligns
  gptm<-round(nr_aligns*gptm_fac,3)
  rpmg<-round(nr_aligns*rpmg_fac,3)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  # Construct return object
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  dt<-data.frame(id=1:n,seqid=seqid,lstart=lstart,lend=lend,rstart=rstart,rend=rend,gaplen=gaplen,
                 strand=strand,nAligns=nr_aligns,gptm=gptm,rpmg=rpmg,nProbes=nProbes)
  ga<-new("gapSites")
  ga@dt<-dt[order(dt$seqid,dt$lend,dt$rstart),]
  ga@nAligns<-round(nAligns,0)
  ga@nAlignGaps<-round(nAlignGaps,0)
  return(ga)
}


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
#  Accessors + Operators for gapSites, dnaGapSites and aaGapSites classes
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Accessors for table column vectors

setMethod("seqid" ,    "gapSites",function(x)      x@dt$seqid)
setMethod("lstart",    "gapSites",function(x)      x@dt$lstart)
setMethod("lend"  ,    "gapSites",function(x)      x@dt$lend)
setMethod("rstart",    "gapSites",function(x)      x@dt$rstart)
setMethod("rend"  ,    "gapSites",function(x)      x@dt$rend)
setMethod("strand",    "gapSites",function(x,...)  x@dt$strand)
setMethod("nAligns",   "gapSites",function(object) object@nAligns)
setMethod("nAlignGaps","gapSites",function(object) object@nAlignGaps)
setMethod("gptm",      "gapSites",function(x)      x@dt$gptm)
setMethod("rpmg",      "gapSites",function(x)      x@dt$rpmg)
setMethod("seqnames",  "gapSites",function(x)      levels(x@dt$seqid))
setMethod("getProfile","gapSites",function(x)      x@profile)

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  "[" Operator

# setMethod("[","gapSites",function(x,i,j,drop){
#   ga<-new("gapSites")
#   ga@dt<-x@dt[i,j]
#   ga@nAligns<-x@nAligns
#   ga@nAlignGaps<-x@nAlignGaps
#   return(ga)
# })

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# getLann, getRann

setMethod("getLann","gapSites",function(x) return(x@lann))
setMethod("getRann","gapSites",function(x) return(x@rann))

setMethod("sortTable","gapSites",function(x){
  o<-order(x@dt$seqid,x@dt$lend,x@dt$rstart)
  cr<-new(class(x))
  cr@dt<-x@dt[o,]
  cr@lann<-x@lann[o,]
  cr@rann<-x@rann[o,]
  cr@nAligns<-x@nAligns
  cr@nAlignGaps<-x@nAlignGaps
  cr@profile<-x@profile
  if(.hasSlot(x,"seq"))
    cr@seq<-x@seq[o]
  return(cr)
})           



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  Merge for gapSites

merge.gapSites<-function(x,y,digits=10,...)
{
  if(!is(y,"gapSites"))
    stop("[merge.gapSites] y must be of class gapSites!")
  message("[merge.gapSites] strand is set to undefined ('*')!")
  
  bm<-Sys.localeconv()[7]
  message("[merge.gapSites] rbind for ",format(x@nAligns,big.mark=bm,digits=digits)
      ," aligns from x and ",format(y@nAligns,big.mark=bm,digits=digits)," aligns from y.")
  
  dt<-rbind(x@dt,y@dt)
  nAligns=as.double(x@nAligns)+as.double(y@nAligns)
  nAlignGaps=as.double(x@nAlignGaps)+as.double(y@nAlignGaps)
  
  sm<-summaryBy(nAligns+nProbes~seqid+lend+rstart,data=dt,FUN=sum,keep.names=TRUE)  
  mx<-summaryBy(rend+gaplen~seqid+lend+rstart,data=dt,FUN=max,keep.names=TRUE)
  mn<-summaryBy(lstart~seqid+lend+rstart,data=dt,FUN=min,keep.names=TRUE)
  
  ga<-gapSites(seqid=sm$seqid,lstart=mn$lstart,
               lend=sm$lend,rstart=sm$rstart,rend=mx$rend,
               gaplen=mx$gaplen,nr_aligns=sm$nAligns,
               nAligns=nAligns,nAlignGaps=nAlignGaps,nProbes=sm$nProbes)
  return(ga)
}

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
# Junction Sites
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("trim_left","gapSites",function(x,maxlen)
{
  if(!is.numeric(maxlen))
    stop("[trim_left.gapRanges] maxlen must be numeric!")
  if(length(maxlen)>1)
    stop("[trim_left.gapRanges] length(maxlen) must be 1!")
  
  maxlen<-as.integer(maxlen)-1L
  cobj<-deparse(substitute(x))
  etxt<-paste(cobj,"@dt$lstart<-pmax(",cobj,"@dt$lstart,",cobj,"@dt$lend-",maxlen,")",sep="")
  eval.parent(parse(text=etxt))
  return(invisible())
})

setMethod("trim_right","gapSites",function(x,maxlen)
{
  if(!is.numeric(maxlen))
    stop("[trim_right.gapRanges] maxlen must be numeric!")
  if(length(maxlen)>1)
    stop("[trim_right.gapRanges] length(maxlen) must be 1!")
  
  maxlen<-as.integer(maxlen)-1L
  cobj<-deparse(substitute(x))
  etxt<-paste(cobj,"@dt$rend<-pmax(",cobj,"@dt$rend,",cobj,"@dt$rstart+",maxlen,")",sep="")
  eval.parent(parse(text=etxt))
  return(invisible())
})

setMethod("resize_left","gapSites",function(x,len)
{
  if(!is.numeric(len))
    stop("[trim_left.gapRanges] len must be numeric!")
  if(length(len)>1)
    stop("[trim_left.gapRanges] length(len) must be 1!")
  len<-as.integer(len)-1L  
  cobj<-deparse(substitute(x))
  etxt<-paste(cobj,"@dt$lstart<-",cobj,"@dt$lend-",len,sep="")
  eval.parent(parse(text=etxt))
  return(invisible())
})

setMethod("resize_right","gapSites",function(x,len)
{
  if(!is.numeric(len))
    stop("[trim_right.gapRanges] len must be numeric!")
  if(length(len)>1)
    stop("[trim_right.gapRanges] length(len) must be 1!")
  len<-as.integer(len)-1L
  cobj<-deparse(substitute(x))
  etxt<-paste(cobj,"@dt$rend<-",cobj,"@dt$rstart+",len,sep="")
  eval.parent(parse(text=etxt))
  return(invisible())
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
# lJunc function
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("lJunc","gapSites",function(x,featlen,gaplen,unique=FALSE,strand,...){
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare featlen and gaplen
  if(!is.numeric(featlen))
    stop("[lJunc.gapRanges] featlen must be numeric!")
  if(length(featlen)>1)
    stop("[lJunc.gapRanges] featlen must have length 1!")
  if(!is.numeric(gaplen))
    stop("[lJunc.gapRanges] gaplen must be numeric!")
  if(length(gaplen)>1)
    stop("[lJunc.gapRanges] gaplen must have length 1!")

  featlen<-as.integer(featlen)
  gaplen<-as.integer(gaplen)
  # Position: 1-based last-exon position
  llen<-featlen-1L
  
  # + + + + + + + + + + + + + + + + + + #
  #  featlen        gaplen
  #        |        |
  #        4321 12345
  #   feat #### ----- gap
  #           |
  #           lend
  #
  # position=featlen
  # + + + + + + + + + + + + + + + + + + #
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare strand
  if(missing(strand))
    stop("[lJunc.gapRanges] strand argument is not optional!")
  if(!is.character(strand))
    stop("[lJunc.gapRanges] strand must be character:",strand)
  if(length(strand)>1)
    stop("[lJunc.gapRanges] strand must have length 1!")
  if(!is.element(strand,strandlevels))
    stop("[lJunc.gapRanges] strand must be valid strand level (+,-,*)!")
  
  cr<-new("cRanges",seqid=x@dt$seqid,start=x@dt$lend-llen,width=featlen+gaplen,strand=strand,position=featlen,id=x@dt$id)
  cr@dt$nAligns<-x@dt$nAligns
  
  # Remove multiplets which arise from alternative (right) sites
  if(unique)
  {
    # Group table
    smin<-summaryBy(id+position~seqid+start+end,data=cr@dt,FUN=min,keep.names=TRUE)    
    ssum<-summaryBy(nAligns~seqid+start+end,data=cr@dt,FUN=sum,keep.names=TRUE)
    slen<-summaryBy(id~seqid+start+end,data=cr@dt,FUN=length,keep.names=TRUE)
    # Copy strand (cannot be grouped)    
    dtu<-cbind(smin[,1:3],data.frame(strand=cr@dt$strand[smin$id],
                    position=smin$position,id=smin$id,nAligns=ssum$nAligns))
    cr@dt<-dtu
  }
  return(cr)
})


setMethod("rJunc","gapSites",function(x,featlen,gaplen,unique=FALSE,strand,...){
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare featlen and gaplen
  if(!is.numeric(featlen))
    stop("[rJunc.gapRanges] featlen must be numeric!")
  if(length(featlen)>1)
    stop("[rJunc.gapRanges] featlen must have length 1!")
  if(!is.numeric(gaplen))
    stop("[rJunc.gapRanges] gaplen must be numeric!")
  if(length(gaplen)>1)
    stop("[rJunc.gapRanges] gaplen must have length 1!")
  featlen<-as.integer(featlen)
  gaplen<-as.integer(gaplen)
  
  # + + + + + + + + + + + + + + + + + + #
  #   gaplen        featlen
  #        |        |
  #        4321 12345
  #   gap  ---- ##### feat
  #             |
  #             rstart
  #
  #
  # + + + + + + + + + + + + + + + + + + #
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare strand
  if(missing(strand))
    stop("[rJunc.gapRanges] strand argument is not optional!")
  if(!is.character(strand))
    stop("[rJunc.gapRanges] strand must be character!")
  if(length(strand)>1)
    stop("[rJunc.gapRanges] strand must have length 1!")
  if(!is.element(strand,strandlevels))
      stop("[rJunc.gapRanges] strand must be valid strand level (+,-,*)!") 
  
  cr<-new("cRanges")
  cr@dt<-data.frame(seqid=x@dt$seqid,start=x@dt$rstart-gaplen,end=x@dt$rstart+featlen-1L,strand=strand,position=featlen,id=x@dt$id)
  cr@dt$nAligns<-x@dt$nAligns
  
  # Remove multiplets which arise from alternative (right) sites
  if(unique)
  {
    smin<-summaryBy(id+position~seqid+start+end,data=cr@dt,FUN=min,keep.names=TRUE)    
    ssum<-summaryBy(nAligns~seqid+start+end,data=cr@dt,FUN=sum,keep.names=TRUE)
    slen<-summaryBy(id~seqid+start+end,data=cr@dt,FUN=length,keep.names=TRUE)
    # Copy strand (cannot be grouped)    
    dtu<-cbind(smin[,1:3],data.frame(strand=cr@dt$strand[smin$id],
                                     position=smin$position,id=smin$id,nAligns=ssum$nAligns))
    cr@dt<-dtu
  }
  return(cr)
})


setMethod("lrJunc","gapSites",function(x,lfeatlen,rfeatlen,strand,...)
{
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare featlen and gaplen
  if(missing(lfeatlen))
    stop("[lrJunc.gapSites] lfeatlen is not optional!")
  if(!is.numeric(lfeatlen))
    stop("[lrJunc.gapSites] lfeatlen must be numeric!")
  if(length(lfeatlen)>1)
    stop("[lrJunc.gapSites] lfeatlen must have length 1!")
  if(missing(rfeatlen))
    stop("[lrJunc.gapSites] rfeatlen is not optional!")
  if(!is.numeric(rfeatlen))
    stop("[lrJunc.gapSites] rfeatlen must be numeric!")
  if(length(rfeatlen)>1)
    stop("[lrJunc.gapSites] rfeatlen must have length 1!")
  lfeatlen<-as.integer(lfeatlen)
  rfeatlen<-as.integer(rfeatlen)

  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare strand
  if(missing(strand))
    stop("[lrJunc.gapSites] strand argument is mandatory!")
  if(length(strand)>1)
    stop("[lrJunc.gapSites] strand must be one unique value!")
  
  ga<-gapSites(seqid=x@dt$seqid,
                lstart=x@dt$lend-lfeatlen+1L,
                lend=x@dt$lend,
                rstart=x@dt$rstart,
                rend=x@dt$rstart+rfeatlen-1L,
                gaplen=x@dt$gaplen,
                strand=strand,
                nr_aligns=x@dt$nAligns,
                nAligns=x@nAligns,
                nAlignGaps=x@nAlignGaps)

  ga@dt$lfeatlen<-lfeatlen
  ga@dt$rfeatlen<-rfeatlen
  return(ga)
})

#setGeneric("lJuncStrand",    function(x,featlen,gaplen,...)standardGeneric("lJuncStrand"))
setMethod("lJuncStrand","gapSites",function(x,featlen,gaplen,...){
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare featlen and gaplen
  if(!is.numeric(featlen))
    stop("[lJunc.gapRanges] featlen must be numeric!")
  if(length(featlen)>1)
    stop("[lJunc.gapRanges] featlen must have length 1!")
  if(!is.numeric(gaplen))
    stop("[lJunc.gapRanges] gaplen must be numeric!")
  if(length(gaplen)>1)
    stop("[lJunc.gapRanges] gaplen must have length 1!")
  
  featlen<-as.integer(featlen)
  gaplen<-as.integer(gaplen)
  # Position: 1-based last-exon position
  llen<-featlen-1L
  
  # + + + + + + + + + + + + + + + + + + #
  #  featlen        gaplen
  #        |        |
  #        4321 12345
  #   feat #### ----- gap
  #           |
  #           lend
  #
  # position=featlen
  # + + + + + + + + + + + + + + + + + + #
  
  cr<-new("cRanges",seqid=x@dt$seqid,start=x@dt$lend-llen,width=featlen+gaplen,strand=x@dt$strand,position=featlen,id=x@dt$id)
  cr@dt$nAligns<-x@dt$nAligns
  return(cr)
})


# setGeneric("rJuncStrand",    function(x,featlen,gaplen,...)standardGeneric("rJuncStrand"))
setMethod("rJuncStrand","gapSites",function(x,featlen,gaplen,...){
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare featlen and gaplen
  if(!is.numeric(featlen))
    stop("[rJunc.gapRanges] featlen must be numeric!")
  if(length(featlen)>1)
    stop("[rJunc.gapRanges] featlen must have length 1!")
  if(!is.numeric(gaplen))
    stop("[rJunc.gapRanges] gaplen must be numeric!")
  if(length(gaplen)>1)
    stop("[rJunc.gapRanges] gaplen must have length 1!")
  featlen<-as.integer(featlen)
  gaplen<-as.integer(gaplen)
  
  # + + + + + + + + + + + + + + + + + + #
  #   gaplen        featlen
  #        |        |
  #        4321 12345
  #   gap  ---- ##### feat
  #             |
  #             rstart
  #
  #
  # + + + + + + + + + + + + + + + + + + #
  
  cr<-new("cRanges")
  cr@dt<-data.frame(seqid=x@dt$seqid,start=x@dt$rstart-gaplen,end=x@dt$rstart+featlen-1L,strand=x@dt$strand,position=featlen,id=x@dt$id)
  cr@dt$nAligns<-x@dt$nAligns
  return(cr)
})


setMethod("lrJuncStrand","gapSites",function(x,lfeatlen,rfeatlen,...)
{
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Prepare featlen and gaplen
  if(missing(lfeatlen))
    stop("[lrJunc.gapSites] lfeatlen is not optional!")
  if(!is.numeric(lfeatlen))
    stop("[lrJunc.gapSites] lfeatlen must be numeric!")
  if(length(lfeatlen)>1)
    stop("[lrJunc.gapSites] lfeatlen must have length 1!")
  if(missing(rfeatlen))
    stop("[lrJunc.gapSites] rfeatlen is not optional!")
  if(!is.numeric(rfeatlen))
    stop("[lrJunc.gapSites] rfeatlen must be numeric!")
  if(length(rfeatlen)>1)
    stop("[lrJunc.gapSites] rfeatlen must have length 1!")
  lfeatlen<-as.integer(lfeatlen)
  rfeatlen<-as.integer(rfeatlen)
  
  ga<-gapSites(seqid=x@dt$seqid,
               lstart=x@dt$lend-lfeatlen+1L,
               lend=x@dt$lend,
               rstart=x@dt$rstart,
               rend=x@dt$rstart+rfeatlen-1L,
               gaplen=x@dt$gaplen,
               strand=x@dt$strand,
               nr_aligns=x@dt$nAligns,
               nAligns=x@nAligns,
               nAlignGaps=x@nAlignGaps)
  
  ga@dt$lfeatlen<-lfeatlen
  ga@dt$rfeatlen<-rfeatlen
  return(ga)
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
#  lrCodons (for gapSites)
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("lrCodons","gapSites",function(x,frame=1L){
  
  if(!is.element(frame,1:3))
    stop("[lrCodons.gapSites] frame must be one of 1,2,3!")
  # Make frame 0-based
  nframe<-as.integer(frame)-1L  
  
  if(sum(table(x@dt$strand)>0)>1)
    stop("[lrCodons.gapSites] Strand must be homogenous. Use 'lrJunc'!")
  
  if(x@dt$strand[1]=='+')
  {
    # We're reading from left to right, so  
    # shift left and truncateSeq right to get full codons in frame
    # Shift left side
    lstart<-x@dt$lstart+nframe
    lfeatlen<-x@dt$lfeatlen-nframe
    
    if(any(lstart>=x@dt$lend))
      stop("[lrCodons.gapSites] Left feature size <=0 in frame",frame,"and (+)-strand!")
    
    # truncateSeq right side
    diff_width<-(lfeatlen+x@dt$rfeatlen)%%3
    rend<-x@dt$rend-diff_width
    rfeatlen<-x@dt$rfeatlen-diff_width
    if(any(x@dt$rstart>=rend))
      stop("[lrCodons.gapSites] Right feature size <=0 after truncation to codon in frame",frame,"and (+)-strand!")
  }
  else
  {
    # We're reading from right to left, so  
    # shift right and truncateSeq left to get full codons in frame
    # Shift right
    
    # truncateSeq left
    rend<-x@dt$rend-nframe   
    rfeatlen<-x@dt$rfeatlen-nframe       
    diff_width<-(rfeatlen+x@dt$lfeatlen)%%3
    
    if(any(rend<=x@dt$rstart))
      stop("[lrCodons.gapSites] Right feature size <=0 in frame",frame,"and (-)-strand!")    
    
    
    lstart<-x@dt$lstart+diff_width
    lfeatlen<-x@dt$lfeatlen-diff_width
    if(any(lstart>=x@dt$lend))
      stop("[lrCodons.gapSites] Left feature size <=0 after truncation to codon in frame",frame,"and (-)-strand!")
  }
  
  dt<-data.frame(
    id=1:length(lstart),
    seqid=x@dt$seqid,
    lstart=lstart,lend=x@dt$lend,
    rstart=x@dt$rstart,rend=rend,
    gaplen=x@dt$gaplen,strand=x@dt$strand,
    nAligns=x@dt$nAligns,
    gptm=x@dt$gptm,rpmg=x@dt$rpmg,
    nProbes=x@dt$nProbes,
    lfeatlen=lfeatlen,rfeatlen=rfeatlen,
    frame=frame
  )
  
  # copy featlen's from input
  ga<-new("gapSites")
  ga@dt<-dt
  ga@nAligns<-x@nAligns
  ga@nAlignGaps<-x@nAlignGaps
  return(ga)
})



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#
#  dnaGapSites : Add dna sequence to gapSites object
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("dnaGapSites","gapSites",function(x,dnaset,strand){
  if(!is(dnaset,"DNAStringSet"))
    stop("[dnaRanges.gapSites] dnaset must be DNAStringSet!\n")
  dnames<-names(dnaset)
  if(is.null(dnames))
    stop("[dnaRanges.gapSites] names(dnaset) must not be NULL!\n")
  
  sqnames<-seqid(x)
  n<-length(dnames)
  seqs<-DNAStringSet()
  lvl<-RleViewsList()
  rvl<-RleViewsList()
  
  # Create empty template
  ans<-x@dt[0]
  j<-1
  for(i in 1:n)
  {
    lg<-match(dnames[i],sqnames)
    if(!is.na(lg))
    {
      dt<-x@dt[x@dt$seqid==dnames[i],]
      # Add data for seqid to data.frame
      ans<-rbind(ans,dt)
      # Add seqname
      #ans_sqnames[j]<-dnames[i]
      # Read dna-ranges
      lvl[[j]]<-Views(dnaset[[i]],start=dt$lstart,end=dt$lend)
      rvl[[j]]<-Views(dnaset[[i]],start=dt$rstart,end=dt$rend)
      j<-j+1
    }
  }
  # chain up dna
  lds<-DNAStringSet(unlist(lapply(lvl,as.character)))
  rds<-DNAStringSet(unlist(lapply(rvl,as.character)))
  dss<-xscat(lds,rds)
  # Remove unused factor levels
  ans$seqid<-factor(ans$seqid)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  strand
  nrows<-nrow(ans)
  if(missing(strand))
  {
    mtc<-match(ans$id,x@dt$id)
    ans$strand<-x@dt$strand[mtc]
    revStrand<-ans$strand=="-"
    dss[revStrand]<-reverseComplement(dss[revStrand])
    gap_pos<-nchar(lds)    
  } else
  {
    if(strand=="-"|| strand==0)
    {
      ans$strand<-factor(rep("-",nrows),levels=strandlevels)
      dss<-reverseComplement(dss)
      gap_pos<-nchar(lds)
    } else
    {
      ans$strand<-factor(rep("+",nrows),levels=strandlevels)
      gap_pos<-nchar(rds)
    }
  }
  
  # Assembly of result
  o<-order(ans$seqid,ans$lend,ans$rstart)
  
  ga<-new("dnaGapSites")
  ga@dt<-ans[o,]
  ga@nAligns<-x@nAligns
  ga@nAlignGaps<-x@nAlignGaps
  # Position = 1-based postion of last exon nucleotide
  ga@dt$position<-gap_pos
  ga@seq<-dss[o]
  return(ga)
})



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#                                                                                                   #
#  Proteomic functions: translate and trypsinCleave                                                 #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("translate","dnaGapSites",function(x){
  
  # + + + + + + + + + + + + + + + + + + #
  #  Exclude data with N's in Sequence
  af<-alphabetFrequency(x@seq)
  withN<-af[,"N"]>0
  noN<-!withN
  
  # + + + + + + + + + + + + + + + + + + #
  # Create new object
  aa<-new("aaGapSites")
  aa@nAligns<-x@nAligns
  aa@nAlignGaps<-x@nAlignGaps 
  
  if(sum(withN)>0)
  {
    message("[translate.dnaGapSites] Excluding ",sum(withN),"/",nrow(x@dt)," rows because of N's!")
    aa@dt<-x@dt[noN,]
    aa@seq<-translate(x@seq[noN])
  }
  else {
    aa@dt<-x@dt
    aa@seq<-translate(x@seq)
  }
    
  # + + + + + + + + + + + + + + + + + + #
  # Range:    1 2 3 | 4 5 6 | 7 8 9
  # Position:         4
  # Corrected position: 2 (4L %/% 3L = 1L)
  # + + + + + + + + + + + + + + + + + + #
  aa@dt$position<-(aa@dt$position %/% 3L)+1L

  return(aa)
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  seqlogo

setMethod("seqlogo","dnaGapSites",function(x,strand="+",useStrand=FALSE,...){
  if(length(x@seq)==0)
    stop("[seqlogo.cdRanges] length(DNAStringSet)=0!")
  # Generates seqLogo for whole DNAStringSet
  if(useStrand)
  {
    strandseqs<-x@dt$strand==strand
    if(sum(strandseqs)==0)
      stop("[seqlogo.cdRanges] No range found for strand '",strand,"'!")
    cs<-consensusMatrix(x@seq[strandseqs],as.prob=T,baseOnly=TRUE)
  }
  else
    cs<-consensusMatrix(x@seq,as.prob=T,baseOnly=TRUE)
  pwm<-makePWM(cs[1:4,])
  seqLogo(pwm,ic.scale=F)
})

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

setMethod("write.files","aaGapSites",function(x,path,filename,...){
  tablefile<-file.path(path,paste(filename,".csv",sep=""))
  fastafile<-file.path(path,paste(filename,".fa",sep=""))
  df<-x@dt
  fasta_id<-1:nrow(df)
  df$fasta_id<-fasta_id
  df$seq<-as.character(x@seq)
  seq<-x@seq
  names(seq)<-paste(fasta_id,filename,sep="|")
  write.table(df,tablefile,sep=";",dec=",",row.names=F,...)
  writeXStringSet(seq,fastafile)
  return(invisible())
})

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#                                                                                                   #
#  Read gapAlings and merged gapSites from BAM via rbamtools                                        #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

# Not exported
do_group_align_data<-function(dt,nAligns,nAlignGaps,startid=1)
{
  # Expects a data.frame with columns:
  # seqid,lstart,lend,rstart,lend,nAligns
  
  # Do grouping by lend,rstart
  dtlen<-summaryBy(position~seqid+lend+rstart,dt,FUN=length,keep.names=TRUE)
  dtmin<-summaryBy(lstart~seqid+lend+rstart,dt,FUN=min,keep.names=TRUE)
  dtmax<-summaryBy(rend~seqid+lend+rstart,dt,FUN=max,keep.names=TRUE)  
  n<-nrow(dtlen)
  # RESET id + strand:
  sites<-data.frame(id=1:n,seqid=dtlen$seqid,lstart=dtmin$lstart,
                    lend=dtlen$lend,rstart=dtlen$rstart,rend=dtmax$rend,
                    gaplen=dtlen$rstart-dtlen$lend-1L,nAligns=dtlen$position,
                    strand=factor(rep("*",n),levels=strandlevels))

  # Calculate new columns for full gapSites object:
  rpmg_fac<-1e6/nAlignGaps
  gptm_fac<-1e7/nAligns
  sites$gaplen<-sites$rstart-sites$lend-1L
  sites$gptm<-round(sites$nAligns*gptm_fac,3)
  sites$rpmg<-round(sites$nAligns*rpmg_fac,3)
  return(invisible(sites))  
}


getGapSites<-function(reader,seqid,startid=1)
{
  if(missing(reader))
    stop("[getGapSites] reader argument is not optional!")
  if(!is(reader,"bamReader"))
    stop("[getGapSites] reader must be bamReader!")
  if(!index.initialized(reader))
    stop("[getGapSites] reader must have initialized Index!")
  if(!is.numeric(seqid))
    stop("[getGapSites] seqid must be numeric!")
  if(length(seqid)!=1)
    stop("[getGapSites] seqid must have length 1")
  if(seqid<=0)
    stop("[getGapSites] seqid must be positive")
  
  ref<-getRefData(reader)
  n<-nrow(ref)
  if(seqid>n)
    stop("[getGapSites] seqid '",seqid,"' must be <=",n,"!")
  
  # Read gaps from BAM
  coords<-c(ref$ID[seqid],0,ref$LN[seqid])
  gl<-gapList(reader,coords)
  
  # BAM-file may contain no reads for given ref
  if(size(gl)==0)
      return(invisible(new("gapSites")))
  
  # Calculate needed columns
  dt<-data.frame(as.data.frame.gapList(gl))
  names(dt)[c(1,5,7)]<- c("seqid","lend","rstart")
  
  dt$lstart<-dt$lend-dt$left_cigar_len+1
  dt$rend<-dt$rstart+dt$right_cigar_len-1
  sites<-do_group_align_data(dt,nAligns(gl),nAlignGaps(gl),startid=startid)
  nRows<-nrow(sites)
  sites$seqid<-factor(rep(ref$SN[seqid],nRows))
    
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  create gapSites Object
  ga<-new("gapSites")
  ga@dt<-sites
  ga@dt$nProbes<-1
  ga@nAligns<-nAligns(gl)
  ga@nAlignGaps<-nAlignGaps(gl)
  return(invisible(ga))
}

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  alignGapProbes: Reads align gaps from entire bam file.
alignGapProbes<-function(reader,startid=1)
{
  if(missing(reader))
    stop("[alignGapProbes] reader argument is not optional!")
  if(!is(reader,"bamReader"))
    stop("[alignGapProbes] reader must be bamReader!")
  if(!rbamtools::isOpen(reader))
    stop("[alignGapProbes] reader must be open!")
  if(!index.initialized(reader))
    stop("[alignGapProbes] reader must have initialized index!")
  l<-list()
  ref<-getRefData(reader)
  n<-nrow(ref)
  
  bm<-Sys.localeconv()[7] 
  nAligns<-0
  k<-1
  for(i in 1:n)
  {
    gal<-getGapSites(reader=reader,seqid=i,startid=startid)
    # 
    if(nrow(gal@dt)>0)
    {
      l[[k]]<-gal
      nAligns<-nAligns+nAligns(l[[k]])
      startid=max(l[[k]]@dt$id)+1
      k<-k+1
      message("[alignGapProbes] seqid: ",format(i,width=2),"\tnAligns:",format(nAligns,big.mark=bm,width=10))
    }
    else
      message("[alignGapProbes] seqid: ",format(i,width=2),"\tnAligns:",format(0,big.mark=bm,width=10))
  }

  # +++++++++++++++++++++++++++++++++++
  ga<-new("gapSites")
  ga@dt<-do.call(rbind,lapply(l,getDataFrame))
  ga@dt$nProbes<-1
  # Push nProbes behind nAligns
  ga@dt<-ga@dt[,c(1:9,12,10,11)]
  ga@nAligns<-sum(unlist(lapply(l,nAligns)))
  ga@nAlignGaps<-sum(unlist(lapply(l,nAlignGaps)))
  
  # +++++++++++++++++++++++++++++++++++
  #  rpmg and gptm
  rpmg_fac<-1e6/ga@nAlignGaps
  gptm_fac<-1e7/ga@nAligns
  ga@dt$gptm<-round(ga@dt$nAligns*gptm_fac,3)
  ga@dt$rpmg<-round(ga@dt$nAligns*rpmg_fac,3)
  return(invisible(ga))
}

alignGapList<-function(reader,digits=3)
{
  if(missing(reader))
    stop("[alignGapList] reader argument is not optional!")
  if(!is(reader,"bamReader"))
    stop("[alignGapList] reader must be bamReader!")
  if(!rbamtools::isOpen(reader))
    stop("[alignGapList] reader must be open!")
  if(!index.initialized(reader))
    stop("[alignGapList] reader must have initialized index!")

  bm<-Sys.localeconv()[7]
  bgl<-bamGapList(reader)
  dtb<-as.data.frame(bgl)
  n<-nrow(dtb)
  dtb$strand<-factor(rep("*",n),levels=strandlevels)
  
  # +++++++++++++++++++++++++++++++++++
  #  rpmg and gptm
  n_aligns<-nAligns(bgl)
  n_gap_aligns<-nAlignGaps(bgl)
  rpmg_fac<-1e6/n_gap_aligns
  gptm_fac<-1e7/n_aligns
  dtb$gptm<-round(dtb$nAligns*gptm_fac,digits=digits)
  dtb$rpmg<-round(dtb$nAligns*rpmg_fac,digits=digits)

  # +++++++++++++++++++++++++++++++++++
  # Create gapSites object
  ga<-new("gapSites")
  ga@nAligns<-n_aligns
  ga@nAlignGaps<-n_gap_aligns
  ga@dt<-dtb[order(dtb$seqid,dtb$lend,dtb$rstart),]
  return(ga)
}

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

readMergedBamGaps<-function(infiles,idxInfiles=paste(infiles,".bai",sep=""),digits=3)
{
  bm<-Sys.localeconv()[7]
  n<-length(infiles)
  
  profile<-data.frame(infile=infiles)
  profile$nAligns<-0      # nAligns for each BAM-file
  profile$nAlignGaps<-0   # nAlignGaps for each BAM-file
  profile$nSites<-0       # number of sites for each BAM-file
  profile$cSites<-0       # number of collected sites (after merging)
    
  for(i in 1:n)
  {
    bam<-infiles[i]
    reader<-bamReader(bam)
    if(!file.exists(idxInfiles[i]))
    {
      message("[readMergedBamGaps] Creating BAM-index.",appendLF=FALSE)
      create.index(reader,idxInfiles[i])
      message("Finished.")
    }
    load.index(reader,idxInfiles[i])
       
    message("[readMergedBamGaps] i:(",format(i,width=2),"/",n,")",appendLF=FALSE)   
    if(i==1)
    {
      ga<-bamGapList(reader)
      profile$nAligns[1]<-nAligns(ga)
      profile$nAlignGaps[1]<-nAlignGaps(ga)
      profile$nSites[1]<-size(ga)
      profile$cSites[1]<-size(ga)
    }
    else
    {
      ga1<-bamGapList(reader)
      ga<-merge.bamGapList(ga,ga1)
      
      profile$nAligns[i]<-nAligns(ga1)
      profile$nAlignGaps[i]<-nAlignGaps(ga1)
      profile$nSites[i]<-size(ga1)
      profile$cSites[i]<-size(ga)
    }
    message("\tnr sites: ",format(size(ga),big.mark=bm,width=9))
  }
  dtb<-as.data.frame(ga)
  n<-nrow(dtb)
  dtb$strand<-factor(rep("*",n),levels=strandlevels)
  
  # +++++++++++++++++++++++++++++++++++
  #  rpmg and gptm
  n_aligns<-nAligns(ga)
  n_gap_aligns<-nAlignGaps(ga)
  rpmg_fac<-1e6/n_gap_aligns
  gptm_fac<-1e7/n_aligns
  dtb$gptm<-round(dtb$nAligns*gptm_fac,digits=digits)
  dtb$rpmg<-round(dtb$nAligns*rpmg_fac,digits=digits)
  
  gap<-new("gapSites")
  # push gqs,nlstart,qmm,nMcs to end of table 
  gap@dt<-dtb[order(dtb$seqid,dtb$lend,dtb$rstart),c(1:7,14,8,9,15,16,10:13)]
  gap@nAligns<-n_aligns
  gap@nAlignGaps<-n_gap_aligns
  gap@profile<-profile
  return(gap)
}


readTabledBamGaps<-function(infiles,idxInfiles=paste(infiles,".bai",sep=""),prof,rpmg=TRUE)
{
  # Check types
  if(!is.character(infiles))
    stop("[readTabledBamGaps] infiles must be character!")
  if(!is.character(idxInfiles))
    stop("[readTabledBamGaps] idxInfiles must be character!")
  if(!is.data.frame(prof))
    stop("[readTabledBamGaps] prof must be data.frame!")
  if(!is.logical(rpmg))
    stop("[readTabledBamGaps] rpmg must be logical!")
  
  # Check sizes
  n<-length(infiles)
  if(nrow(prof)!=n)
    stop("[readTabledBamGaps] Length of infiles must be equal to row number in prof!")
  if(length(idxInfiles)!=n)
    stop("[readTabledBamGaps] Length of infiles must be equal length of idxInfiles!")
  
  if(!all(lapply(prof,class)=="factor"))
    stop("[readTabledBamGaps] All columns of 'prof' (profile) must be factor!")
  
  # Remove unused factor levels
  for(i in 1:ncol(prof))
    prof[,i]<-factor(prof[,i])
  
  # Create profile table
  profile<-prof
  
  profile$nAligns<-0      # nAligns for each BAM-file
  profile$nAlignGaps<-0   # nAlignGaps for each BAM-file
  profile$nSites<-0       # number of sites for each BAM-file
  profile$cSites<-0       # number of collected sites (after merging)
  
  profile$infile<-infiles
  
  bm<-Sys.localeconv()[7]
  for(i in 1:n)
  {
    bam<-infiles[i]
    reader<-bamReader(bam)
    if(!file.exists(idxInfiles[i]))
    {
      message("[readTabledBamGaps] Creating BAM-index.",appendLF=FALSE)
      create.index(reader,idxInfiles[i])
      message("Finished.")
    }
    load.index(reader,idxInfiles[i])
    message("[readTabledBamGaps] i:(",format(i,width=2),"/",n,")",appendLF=FALSE)    
    if(i==1)
    {
      ga<-bamGapList(reader)
      dt<-as.data.frame(ga)
      kpc<-new("keyProfiler",keyTable=dt[,c("seqid","lend","rstart")],prof=prof,index=1,unique=TRUE)
      kpa<-new("keyProfiler",keyTable=dt[,c("seqid","lend","rstart")],prof=prof,index=1,values=dt$nAligns,unique=TRUE)
      profile$nAligns[1]<-nAligns(ga)
      profile$nAlignGaps[1]<-nAlignGaps(ga)
      profile$nSites[1]<-size(ga)
      profile$cSites[1]<-size(ga)
    }
    else
    {
      ga1<-bamGapList(reader)
      dt<-as.data.frame(ga1)
      addKeyTable(kpc,keyTable=dt[,c("seqid","lend","rstart")],index=i)
      addKeyTable(kpa,keyTable=dt[,c("seqid","lend","rstart")],index=i,values=dt$nAligns)
      profile$nAligns[i]<-nAligns(ga1)
      profile$nAlignGaps[i]<-nAlignGaps(ga1)
      profile$nSites[i]<-size(ga1)
      
      # Merge
      ga<-merge.bamGapList(ga,ga1)
      profile$cSites[i]<-size(ga)      
    }
    message("\tnr sites: ",format(size(ga),big.mark=bm,width=9))
  }
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Create data.frame from bamGapList
  dtb<-as.data.frame(ga)
  n<-nrow(dtb)
  dtb$strand<-factor(rep("*",n),levels=strandlevels)
  
  # +++++++++++++++++++++++++++++++++++
  #  rpmg and gptm
  n_aligns<-nAligns(ga)
  n_gap_aligns<-nAlignGaps(ga)
  rpmg_fac<-1e6/n_gap_aligns
  gptm_fac<-1e7/n_aligns
  dtb$gptm<-round(dtb$nAligns*gptm_fac,3)
  dtb$rpmg<-round(dtb$nAligns*rpmg_fac,3)
  #setkeyv(dtb,c("seqid","lend","rstart"))
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Append KeyTable values
  res<-appendKeyTable(kpc,dtb,prefix="c.")
  res<-appendKeyTable(kpa,res,prefix="aln.")
  
  if(rpmg)
    res<-appendKeyTable(kpa,res,prefix="rpmg.",rateFactor=1e6,digits=3)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Create gapSites object
  gap<-new("gapSites")
  gap@nAligns<-n_aligns
  gap@nAlignGaps<-n_gap_aligns
  gap@profile<-profile
  
  # Reorder first 16 columns
  nCol<-ncol(res)
  gap@dt<-res[order(res$seqid,res$lend,res$rstart),c(4,1,5,2,3,6,7,14,8:13,15,16,17:nCol)]
  return(gap)
}



rangeByGeneName<-function(reader,genome,gene,complex=TRUE)
{
  if(!is(reader,"bamReader"))
    stop("[rangeByGeneName] 'reader' must be of class 'bamReader'!")
  if(!is(genome,"refGenome"))
     stop("[rangeByGeneName] 'genome' must be of class 'refGenome'!")
  if(!rbamtools::isOpen(reader))
       stop("[rangeByGeneName] reader must be opened!")     
  if(!index.initialized(reader))
       stop("[rangeByGeneName] reader must have initialized index!")
  if(!is.character(gene))
       stop("[rangeByGeneName] gene must be character!")
  if(length(gene)>1)
       stop("[rangeByGeneName] gene must have length 1!")
  if(!is.logical(complex))
       stop("[rangeByGeneName] complex must be logical!")
  if(length(complex)>1)
       stop("[rangeByGeneName] complex must have length 1!")
     
  exons<-extractByGeneName(genome,gene)
  if(is.null(exons))
  {
     message("[rangeByGeneName] Gene '",gene,"' does not match any annotation. Returning NULL.")
     return(invisible(NULL))
  }
  pos<-getGenePositions(exons)
  seq<-as.character(pos$seqid) 
  start<-as.numeric(pos$start)
  end<-as.numeric(pos$end)
  
  ref<-getRefData(reader)
  mtc<-match(seq,ref$SN)
  if(is.na(mtc))
    stop("[rangeByGeneName] Missing seqid for reference sequence '",seq,"'!")
  seqid<-ref$ID[mtc]
  
  return(bamRange(reader,coords=c(seqid,start,end),complex=complex))
}

countByGeneName<-function(object,infiles,idxInfiles=paste(infiles,".bai",sep=""),gene,tag="N")
{
  if(!is(object,"refGenome"))
    stop("[countByGeneName] 'object' must be of class 'refGenome'!")
  if(missing(gene))
    stop("[countByGeneName] 'gene' argument is not optional!")
  if(length(gene)>1)
    stop("[countByGeneName] Only single genes!")
  if(length(tag)>1)
    stop("[countByGeneName] Only single tags!")
  
  exons<-extractByGeneName(object,gene)
  # Gene does not match?
  if(is.null(exons))
  {
    message("[countByGeneName] Gene '",gene,"' does not match any annotation. Returning NULL.\n")
    return(invisible(NULL))
  }
  
  # Calls 'exon'-specific version of getGenePositions
  # Different for ensemblGenome and ucscGenome (see: 'by' argument)
  pos<-getGenePositions(exons)
  
  # Try to deal with multiple results
  seq<-as.character(unique(pos$seqid))
  if(length(seq)>1)
    stop("[countByGeneName] Found gene on multiple seqid's!")
  
  start<-min(as.numeric(pos$start))
  end<-max(as.numeric(pos$end))
  n<-length(infiles)
  res<-numeric(n)
  
  tag_match<-0
  
  for(i in 1:n)
  {
    message("[countByGeneName] i:(",format(i,width=2),"/",n,")")
    bam<-infiles[i]
    reader<-bamReader(bam)
    ref<-getRefData(reader)
    mtc<-match(seq,ref$SN)
    if(is.na(mtc))
      stop("[countByGeneName] Annotation seqid '",seq,"' not found in BAM-file (",i,")!")
    
    if(!file.exists(idxInfiles[i]))
    {
      message("[countByGeneName] Creating BAM-index.",appendLF=FALSE)
      create.index(reader,idxInfiles[i])
      message("Finished.\n")
    }
    load.index(reader,idxInfiles[i])
    
    count<-bamCount(reader,coords=c(ref$ID[mtc],start,end))
    if(tag_match==0)
    {
      tag_match<-match(tag,names(count))  
      if(is.na(tag_match))
        stop("[countByGeneName] Tag '",tag,"' not found in count values!")      
    }
    res[i]<-count[tag_match]
  }
  message("[countByGeneName] Finished.")
  return(res)
}


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#                                                                                                   #
# Beginn annotation routines                                                                        #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  Overlap and merge gapSites data with annotation (refGenome)

overlap_genome<-function(qry,ref,verbose=FALSE)
{
  # expects qry and ref as data.frames with columns
  # id, start, end, seqid
  #qryRefs<-table(qry$seqid)
  qRefNames<-levels(qry$seqid)
  nQryRefs<-length(qRefNames)
  
  l<-list()
  k<-1
  bm<-Sys.localeconv()[7]
  for(i in 1:nQryRefs)
  {
    qrs<-qry[qry$seqid==qRefNames[i],]
    qrs<-qrs[order(qrs$start),]
    rfs<-ref[ref$seqid==qRefNames[i],]
    rfs<-rfs[order(rfs$start),]
    nq<-nrow(qrs)
    nr<-nrow(rfs)
    if(verbose)
    {
      message("[overlap_genome] i: ",format(i,width=2),"\tseqid: ",format(qRefNames[i],width=4),"\tquery set: ",format(nq,width=7,big.mark=bm),
          "\tref set:",format(nr,width=6,big.mark=bm))
      if(nq==0|nr==0)
        message("\t\tskipped!")
      else
      {
        message("")
        l[[k]]<-overlap(qry=qrs,ref=rfs)
        k<-k+1
      }      
    }
    else # verbose=FALSE
    {
      if(nq>0&nr>0)
      {
        l[[k]]<-overlap(qry=qrs,ref=rfs)
        k<-k+1
      } 
    }
  }
  return(do.call(rbind,l))
}


setMethod("annotate","gapSites",function(object,genome)
{
  extract<-c("id","start","end","seqid")
  if(is(genome,"refGenome"))
  {
    ref<-genome@ev$gtf[,extract]
    reftbl<-genome@ev$gtf
  }
  else if(is(genome,"data.frame"))
  {
    if(any(!is.element(extract,names(genome))))
      stop("[annotate.gapSites] genome table must contain columns 'id','start','end','seqid'")
    ref<-genome[,extract]
    reftbl<-genome
  }
  else
    stop("[annotate.gapSites] genome must be 'refGenome' or 'data.frame'!")
  
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  ## Preparation of reference
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

  mtc<-match(levels(ref$seqid),levels(object@dt$seqid))
  if(any(is.na(mtc)))
    message("[annotate.gapSites] Skipping ",sum(is.na(mtc))," seqid's in genome!") 
  if(sum(!is.na(mtc))==0)
    stop("[annotate.gapSites] No match between gapSites seqid's and genome seqid's (wrong genome?)!")
  l<-list()
    
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  ##  Annotation of left side (lstart, lend)
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  message("[annotate.gapSites] Annotating left side.")
  qry<-object@dt[order(object@dt$seqid,object@dt$lstart),c("id","lstart","lend","seqid")]
  names(qry)[2:3]<-c("start","end")
  
  l$left<-overlap_genome(qry,ref)
  l$left<-merge(l$left,reftbl,by.x="refid",by.y="id",all.x=TRUE)
  names(l$left)<-paste("left",names(l$left),sep="_")
  names(l$left)[names(l$left)=="left_queryid"]<-"id"
  l$left<-l$left[order(l$left$id),]
  
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  ##  Annotation of right side (rstart,rend)
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  message("[annotate.gapSites] Annotating right side.")
  qry<-object@dt[order(object@dt$seqid,object@dt$rstart),c("id","rstart","rend","seqid")]
  names(qry)[2:3]<-c("start","end")
  
  l$right<-overlap_genome(qry,ref)
  l$right<-merge(l$right,reftbl,by.x="refid",by.y="id",all.x=TRUE)
  names(l$right)<-paste("right",names(l$right),sep="_")
  names(l$right)[names(l$right)=="right_queryid"]<-"id"
  l$right<-l$right[order(l$right$id),]
  
  class(l)<-c("list","annotationResult")  
  message("[annotate.gapSites] Finished.")
  return(invisible(l))
})

setMethod("annotate","ExpressionSet",function(object,genome)
{
  # Ensembl Genome? (enpa)
  fdo<-featureData(object)
  fd<-fdo@data
  fd$id<-1:nrow(fd)
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  message("[annotate.ExpressionSet] Annotating left side.")
  # Annotation of left side (lstart, lend)
  qry<-fd[,c("id","lstart","lend","seqid")]
  names(qry)<-c("id","start","end","seqid")
  
  # no need to sort
  ref<-genome@ev$gtf[genome@ev$gtf$feature=="exon",c("id","start","end","seqid")]
    
  ovg<-overlap_genome(qry,ref)
  if(is.null(ovg))
    stop("[annotate.ExpressionSet] No overlap found for left side (wrong reference genome version?)")
  
  lo<-merge(ovg[,c("queryid","refid")],genome@ev$gtf[,c("id","gene_id","transcript_id","gene_name")],by.x="queryid",by.y="id",all.x=TRUE)
  names(lo)<-paste("left",names(lo),sep="_")
  res<-merge(fd,lo,by.x="id",by.y="left_queryid")
    
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Annotation of right side (rstart, rend)
  message("[annotate.ExpressionSet] Annotating right side.")
  
  qry<-res[,c("id","rstart","rend","seqid")]
  names(qry)[2:3]<-c("start","end")
  
  ovg<-overlap_genome(qry,ref)
  if(is.null(ovg))
    stop("[annotate.ExpressionSet] No overlap found for right side (wrong reference genome version?)")
  
  ro<-merge(ovg[,c("queryid","refid")],genome@ev$gtf[,c("id","gene_id","transcript_id","gene_name")],by.x="refid",by.y="id",all.x=TRUE)
  names(ro)<-paste("right",names(ro),sep="_")
  res<-merge(res,ro,by.x="id",by.y="right_queryid")
  
  
  meta<-data.frame(labelDescription=c("identifier","sequence identifier","left end","right start","left start","right end",
                                      "reference table id","left gene id", "left transcript id", "left gene name",
                                      "reference table id","right gene id", "right transcript id", "right gene name"
                                      ),
                   row.names=colnames(res))
  adf<-new("AnnotatedDataFrame",data=res,varMetadata=meta)
  message("[annotate.ExpressionSet] Finished.")
  return(adf)
})


setReplaceMethod("annotation","gapSites",function(object,value)
{
  if(!is(value,"annotationResult"))
    stop("[annotation<-.gapSites] value must be of type 'annotationResult'. Use 'annotate'!")
  object@lann<-value$left
  object@rann<-value$right
  return(object)
})

setMethod("annotation","gapSites",function(object)
{
  if(is.null(object@lann))
    stop("[annotation.gapSites] Annotation table is NULL! Use 'annotation<-' and 'annotate'!")
  
  #  merge annotations with main table (dt)
  dt<-merge(object@dt,object@lann,by="id",all.x=TRUE)
  dt<-merge(dt,object@rann,by="id",all.x=TRUE)
  return(invisible(new("annGapSites",dt=dt,nAligns=object@nAligns,nAlignGaps=object@nAlignGaps)))
})

#setGeneric("getAnnStrand",function(object)standardGeneric("getAnnStrand"))
setMethod("getAnnStrand","gapSites",function(object)
{
  if(is.null(object@lann) || is.null(object@rann))
    stop("[getAnnStrand.gapSites] lann and rann table must not be NULL! Use 'annotate' and 'annotation<-'!")
  # Create data.frame with left annotation strand entries
  leftstrand<-object@lann[,c("id","left_strand")]
  names(leftstrand)[2]<-"lstr"
  levels(leftstrand$lstr)<-c(levels(leftstrand$lstr),"*")
  leftstrand$lstr[is.na(leftstrand$lstr)]<-"*"
  
  # Create data.frame with right annotation strand entries
  rightstrand<-object@rann[,c("id","right_strand")]
  names(rightstrand)[2]<-"rstr"
  levels(rightstrand$rstr)<-c(levels(rightstrand$rstr),"*")
  rightstrand$rstr[is.na(rightstrand$rstr)]<-"*"
  # Merge left and right strand info into one data.frame
  dstr<-merge(leftstrand,rightstrand,by="id")
  
  # Create combined strand information 
  n<-nrow(dstr)
  dstr$strand<-factor(rep("*",n),levels=strandlevels)
  dstr$strand[dstr$lstr=="+" & dstr$rstr=="+"]<-"+"
  dstr$strand[dstr$lstr=="-" & dstr$rstr=="-"]<-"-"
  
  # Remove left-strand and right-strand columns
  dstr$lstr<-NULL
  dstr$rstr<-NULL
  class(dstr)<-c("data.frame","AnnStrandInfo")
  return(invisible(dstr))
})

setReplaceMethod("strand","gapSites",function(x,value)
{
  if(!is(value,"AnnStrandInfo"))
    stop("[strand<-.gapSites] value must be of class 'AnnStrandInfo'. Use 'getAnnStrand'!")
  
  # Overwrite old strand information
  n<-nrow(x@dt)
  x@dt$strand<-factor(rep("*",n),levels=strandlevels)  
  mtc<-match(x@dt$id,value$id)
  nmtc<-!is.na(mtc)
  x@dt$strand[nmtc]<-value$strand[mtc[nmtc]]
  return(x)
})

#setGeneric("addGeneAlignPart",function(x)standardGeneric("addGeneAlignPart"))
setMethod("addGeneAlignPart","gapSites",function(x)
{
  if(is.null(x@lann) | is.null(x@rann))
    stop("[addGeneAlignPart] Annotation missing. Use 'annotation(x)<-annotate(x,...)'!")
  
  bm<-Sys.localeconv()[7]
  
  atb<-merge(x@dt,x@lann[,c("id","left_gene_id","left_gene_name")],by="id")
  atb<-merge(atb,x@rann[,c("id","right_gene_id","right_gene_name")],by="id")
  
  isna<-is.na(atb$left_gene_id)|is.na(atb$right_gene_id)
  message("[addGeneAlignPart] Skipping ",format(sum(isna),big.mark=bm,width=8)," sites because of missing Gene id.")
  atb<-atb[!isna,]
  
  eq<-atb$left_gene_id==atb$right_gene_id
  message("[addGeneAlignPart] Skipping ",format(sum(!eq),big.mark=bm,width=8)," sites because of unequal Gene id.")
  atb<-atb[eq,c("id","nAligns","left_gene_id")]
  
  sal<-tapply(atb$nAligns,atb$left_gene_id,sum)
  tal<-data.frame(gene_aligns=sal,left_gene_id=names(sal))
  
  mtc<-match(atb$left_gene_id,tal$left_gene_id)
  atb$gene_aligns<-tal$gene_aligns[mtc]
  atb$align_part<-round(atb$nAligns/atb$gene_aligns,4)
  
  # Create copy of incoming object
  ga<-new("gapSites")
  ga@dt<-x@dt
  ga@lann<-x@lann
  ga@rann<-x@rann
  ga@nAligns<-x@nAligns
  ga@nAlignGaps<-x@nAlignGaps  
  
  # Add gene_aligns and align_part column
  mtc<-match(x@dt$id,atb$id)
  ga@dt$gene_aligns<-atb$gene_aligns[mtc]
  ga@dt$align_part<-atb$align_part[mtc]

  return(ga)
})

# Plots tabled distance between inner gap-site and annotation borders
setMethod("plot_diff","annGapSites",function(x,n=20)
{
  tbl1l<-table(abs(x@dt$lend-x@dt$left_end))
  tbl1r<-table(abs(x@dt$rstart-x@dt$right_start))
  tbl2l<-table(abs(x@dt$left_rightDiff))
  tbl2r<-table(abs(x@dt$right_leftDiff))

  op<-par(mfrow=c(2,2))
  barplot(tbl1l[1:n],main="Annotation left-gap_site distance")  
  barplot(tbl1r[1:n],main="Annotation right-gap site distance")
  barplot(tbl2l[1:n],main="Left rightDiff")
  barplot(tbl2r[1:n],main="Right leftDiff")
  par(op)
  return(invisible())
})

# End annotation routines
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#                                                                                                   #
# Beginn alternative sites routines                                                                 #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Basic alt_left_ranks and alt_right_ranks functions
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

  ###################################################################
  # Model:                                                          #
  #             lstart     lend            rstart     rend          #
  #             ###############            ###############          #
  #             |  left exon  |   intron   |  right exon |          #
  #                                                                 #
  # alt_left:               alt            group                    #
  #                                                                 #
  ###################################################################

setMethod("alt_left_ranks","gapSites",function(x){  
  # Function core is generic for left and right:
  # alt_left: grouping by rstart   and   alt(ernatives) = lend
  # Extraction of alt-table
  gt<-as.data.frame(x@dt)[,c(1,2,4,5,7)]
  gt<-data.frame(lapply(gt,as.integer))
  gt<-gt[order(gt$seqid,gt$rstart,gt$gaplen),]  
  return(invisible(.Call("alt_group",gt$id,gt$seqid,gt$rstart,gt$gaplen)))
})


  ###################################################################
  # Model:                                                          #
  #             lstart     lend            rstart     rend          #
  #     5' ->   ###############            ###############   -> 3'  #
  #             |  left exon  |   intron   |  right exon |          #
  #                                                                 #
  # alt_right:            group            alt                      #
  #                                                                 #
  ###################################################################

setMethod("alt_right_ranks","gapSites",function(x){  
  # Function core is generic for left and right:
  # alt_left: grouping by rstart   and   alt(ernatives) = lend
  # Extraction of alt-table
  gt<-as.data.frame(x@dt)[,c(1,2,4,5,7)]
  gt<-data.frame(lapply(gt,as.integer))
  gt<-gt[order(gt$seqid,gt$lend,gt$gaplen),]
  return(invisible(.Call("alt_group",gt$id,gt$seqid,gt$lend,gt$gaplen)))
})

setMethod("alt_ranks","gapSites",function(object)
{
  alt_left<-alt_left_ranks(object)
  alt_right<-alt_right_ranks(object)
  res<-merge(alt_left,alt_right,by="id",suffixes=c("_left","_right"))
  res$seqid<-object@dt$seqid
  res$lend<-object@dt$lend
  res$rstart<-object@dt$rstart  
  return(res)
})

setMethod("plot_diff_ranks","gapSites",function(x){
  alt_left<-alt_left_ranks(x)
  alt_right<-alt_right_ranks(x)
  
  tbll<-table(alt_left$gap_diff)  
  tblr<-table(alt_right$gap_diff)
  op<-par(mfrow=c(1,2))
  barplot(tbll[2:50],main="table left_gap_diff")  
  barplot(tblr[2:50],main="table right_gap_diff")  
  par(op)
  return(invisible())
})

# End alternative sites routines
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #




# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  Coerce
as.data.frame.cRanges<-function(x, row.names=NULL, optional=FALSE, ...){return(x@dt)}
as.data.frame.gapSites<-function(x, row.names=NULL, optional=FALSE, ...){
  res<-as.data.frame(x@dt)
  if(!is.null(x@lann))
    res<-merge.data.frame(res,x@lann[,c("id","left_gene_id","left_gene_name")],by="id")
  if(!is.null(x@rann))
    res<-merge.data.frame(res,x@rann[,c("id","right_gene_id","right_gene_name")],by="id")
  return(res)
}

# Not exported:
getDataFrame<-function(x) {return(x@dt)}
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  Export functions: write csv files

#setGeneric("write.annDNA.tables",function(object,dnaset,filename,featlen=3,gaplen=8,sep=";",dec=".",row.names=FALSE)standardGeneric("write.annDNA.tables"))
setMethod("write.annDNA.tables","gapSites",function(object,dnaset,filename,featlen=3,gaplen=8,sep=";",dec=",",row.names=FALSE)
{
  if(!is(dnaset,"DNAStringSet"))
    stop("[write.AnnotatedDNA.tables.gapSites] dnaset must be 'DNAStringSet'!")
  if(is.null(object@lann))
    stop("[write.AnnotatedDNA.tables.gapSites] No annotation table available! Use 'annotation<-'!")
  if(!is.character(filename))
    stop("[write.AnnotatedDNA.tables.gapSites] filename must be 'character'!")
  
  dt<-annotation(object)@dt
  dt<-dt[order(dt$id),]
  
  ljseq<-dnaRanges(lJuncStrand(object,featlen=featlen,gaplen=gaplen),dnaset,removeUnknownStrand=FALSE,verbose=FALSE)
  lsq<-ljseq@dt
  lsq$leftseq<-as.character(ljseq@seq)
  names(lsq)[names(lsq)=="position"]<-"left_position"
  lsq$seqid<-NULL
  lsq$start<-NULL
  lsq$end<-NULL
  lsq$strand<-NULL
  lsq$nAligns<-NULL
  dt<-merge(dt,lsq,by="id",all.x=TRUE)
  
  rjseq<-dnaRanges(rJuncStrand(object,featlen=featlen,gaplen=gaplen),dnaset,removeUnknownStrand=FALSE,verbose=FALSE)
  rsq<-rjseq@dt
  rsq$rightseq<-as.character(rjseq@seq)
  names(lsq)[names(lsq)=="position"]<-"right_position"
  rsq$seqid<-NULL
  rsq$start<-NULL
  rsq$end<-NULL
  rsq$strand<-NULL
  rsq$nAligns<-NULL
  dt<-merge(dt,rsq,by="id",all.x=TRUE)
  
  old_sci<-options()$scipen
  options(scipen=999)
  write.table(dt,file=filename,sep=sep,dec=dec,row.names=row.names,na="")
  options(scipen=old_sci)
  message("[write.AnnotatedDNA.tables.gapSites] Finished.")
  return(invisible())
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  keyProfiler:                                                                                     #
#  counts occurrence of profile (prof) items defined by factor levels in profile table              #
#  in successively added tables by keys defined in object@ev$dt                                     #
#  unique: allows adding keys only once for each profile line                                       #
#  values: when given, conditions are not simply counted but the given values                       #
#          are summed up for each condition.                                                        #
#          values must either be absent or present for all accumulated data                         #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

# setClass("keyProfiler",representation(ev="environment",unique="logical",counted="logical",useValues="logical"))
setMethod("initialize","keyProfiler",
          function(.Object,
                   keyTable,
                   prof,
                   values,
                   unique=TRUE,
                   index=1L)
          {
            .Object@ev=new.env()
            if(!is(prof,"data.frame"))
              stop("[initialize.keyProfiler] prof must be data.frame!")
            n<-nrow(prof)
            if(n<2)
              stop("[initialize.keyProfiler] prof must at least have two rows (=categories), otherwise there's nothing to count!")
            
            if(!all(unlist(lapply(prof,"class"))=="factor"))
              stop("[initialize.keyProfiler] All columns in prof must be of class 'factor'!")
            if(!is.logical(unique))
              stop("[initialize.keyProfiler] unique must be 'logical'!")
            if(length(unique)>1)
              stop("[initialize.keyProfiler] unique must have length 1!")
            if(!is.numeric(index))
              stop("[initialize.keyProfiler] index must be numeric!")
            if(!is.integer(index))
              index<-as.integer(index)
            if(length(index)>1)
              stop("[initialize.keyProfiler] index must have length 1!")
            if(index<0 || index>n)
              stop("[initialize.keyProfiler] index<0 or index>",n," (points to row in prof)!")
            useValues<-!missing(values)
            if(useValues)
            {
              if(!is.numeric(values))
                stop("[initialize.keyProfiler] Given values must be numeric (they are summed up)!")
              if(length(values)!=nrow(keyTable))
                stop("[initialize.keyProfiler] length(values) unequal row-number in keyTable!")
            }
            
            # Reshape profile prof-table for counting
            .Object@ev$prof<-prof
            nProfCols<-ncol(.Object@ev$prof)
            profLevels<-lapply(.Object@ev$prof,levels)
            for(i in 1:nProfCols)
            {
              tag<-names(profLevels)[i]
              colNames<-paste(tag,profLevels[[i]],sep=".")
              nCol<-length(colNames)
              for(j in 1:nCol)
              {
                lgl<-ifelse(.Object@ev$prof[,tag]==profLevels[[i]][j],1,0)
                .Object@ev$prof[,colNames[j]]<-lgl
              }
              .Object@ev$prof[,tag]<-NULL
            }
            
            # Add profiling columns to keytable
            nProfCols<-ncol(.Object@ev$prof)
            colNames<-names(.Object@ev$prof)
            
            kt<-data.frame(keyTable)
            if(useValues)
            {
              # Sum up values
              for(i in 1:nProfCols)
                if(as.integer(.Object@ev$prof[index,i])==1)
                  kt[,colNames[i]]<-values
              else
                kt[,colNames[i]]<-0
            } else
            {
              # Do simple counting
              for(i in 1:nProfCols)
                kt[,colNames[i]]<-as.integer(.Object@ev$prof[index,i])
            }
            
            .Object@ev$dtb<-kt
            .Object@ev$keyTableNames<-names(keyTable)
            .Object@unique<-unique
            .Object@useValues<-useValues
            .Object@ev$groupExpr<-parse(text=paste("dtb[,lapply(.SD,sum),by=list(",paste(names(keyTable),collapse=","),")]",sep=""))
            if(unique)
            {
              .Object@counted<-rep(FALSE,nrow(prof))
              .Object@counted[index]<-TRUE
            }
            return(.Object)
          })


setMethod("addKeyTable","keyProfiler",function(object,keyTable,index,values)
{
  if(object@unique)
  {
    if(object@counted[index])
      stop("[addKeyTable.keyProfiler] For each index, only one table can be added (or set unique=FALSE)!")
    object@counted[index]<-TRUE
  }
  # Add profiling columns to keytable
  nProfCols<-ncol(object@ev$prof)
  colNames<-names(object@ev$prof)
  
  if(object@useValues)
  {
    # Sum up values
    for(i in 1:nProfCols)
      if(as.integer(object@ev$prof[index,i])==1)
        keyTable[,colNames[i]]<-values
    else
      keyTable[,colNames[i]]<-0
  } else
  {
    # Do simple counting
    for(i in 1:nProfCols)
      keyTable[,colNames[i]]<-as.integer(object@ev$prof[index,i])
  }
  dtb<-data.frame(rbind(object@ev$dtb,keyTable))
  object@ev$dtb<-summaryBy(.~seqid+lend+rstart,dtb,FUN=sum,keep.names=TRUE)
  return(invisible())
})

setMethod("getKeyTable","keyProfiler",function(object) {return(object@ev$dtb)})

setMethod("appendKeyTable","keyProfiler",function(object,keytable,prefix,valFactor,rateFactor,digits)
{
  
  if(!missing(valFactor))
  {
    if(!is.numeric(valFactor))
      stop("[appendKeyTable.keyProfiler] valFactor must be numeric!")
    if(length(valFactor)>1)
      stop("[appendKeyTable.keyProfiler] valFactor must have length 1!")
    if(!missing(rateFactor))
      message("[appendKeyTable.keyProfiler] rateFactor omitted because valFactor is present!")
    
    # When a value Factor is given, each profile column
    # will be multiplied with the given value
    keylen<-length(object@ev$keyTableNames)
    tbl_width<-ncol(object@ev$dtb)
    
    if(missing(digits))
    {
      # Walk through profile columns
      # ToDo: Check for valid table size? (keylen+1)<tbl_width ?
      for(i in (keylen+1):tbl_width)
        object@ev$dtb[,i]<-object@ev$dtb[,i]*valFactor       
    }
    else
    {
      if(!is.numeric(digits))
        stop("[appendKeyTable.keyProfiler] digits must be numeric!")
      if(length(digits)>1)
        stop("[appendKeyTable.keyProfiler] digits must have length 1!")
      
      # Walk through profile columns
      for(i in (keylen+1):tbl_width)
        object@ev$dtb[,i]<-round(object@ev$dtb[,i]*valFactor,digits=as.integer(digits))
    }
    
  }
  else if(!missing(rateFactor))
  {
    if(!is.numeric(rateFactor))
      stop("[appendKeyTable.keyProfiler] rateFactor must be numeric!")
    if(length(rateFactor)>1)
      stop("[appendKeyTable.keyProfiler] rateFactor must have length 1!")
    
    # When a rate Factor is given, each profile column will be converted into a rate
    # i.e. each profile value will be divided by its column-sum.
    # In order to get a sensible rate, the profile columns are then multiplied
    # with the rateFactor.
    keylen<-length(object@ev$keyTableNames)
    tbl_width<-ncol(object@ev$dtb)
    
    if(missing(digits))
    {
      # Walk through profile columns
      for(i in (keylen+1):tbl_width)
      {
        col_sum<-sum(object@ev$dtb[,i])
        if(col_sum>0)
        {
          col_factor<-rateFactor/col_sum
          object@ev$dtb[,i]<-object@ev$dtb[,i]*col_factor          
        }
      }       
    }
    else
    {
      if(!is.numeric(digits))
        stop("[appendKeyTable.keyProfiler] digits must be numeric!")
      if(length(digits)>1)
        stop("[appendKeyTable.keyProfiler] digits must have length 1!")
      
      # Walk through profile columns
      for(i in (keylen+1):tbl_width)
      {
        col_sum<-sum(object@ev$dtb[,i])
        if(col_sum>0)
        {
          col_factor<-rateFactor/col_sum
          object@ev$dtb[,i]<-round(object@ev$dtb[,i]*col_factor,digits=as.integer(digits))
        }
      }       
    }
    
  }
  ans<-merge(keytable,object@ev$dtb,by=object@ev$keyTableNames)
  if(!missing(prefix))
  {
    colNames<-names(object@ev$prof)
    mtc<-match(colNames,names(ans))
    names(ans)[mtc]<-paste(prefix,colNames,sep="")
  }
  return(ans)
})

# End keyProfiler
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Read expression sets

padd<-function(text,len=5L,char="#",left=TRUE)
{
  if(!is.character(text))
    stop("'text' must be character!")
  if(!is.integer(len))
    stop("'len' must be integer!")
  if(length(len)>1)
    stop("'len' must have length 1!")
  if(!is.character(char))
    stop("'char' must be character")
  if(!is.logical(left))
    stop("'left' must be logical")
  
  slen<-nchar(text)
  if(any(slen)>len)
    stop("'len' must be >= nchar(text)!")
  
  padstring<-paste(rep(char,len),collapse="")
  padlen<-len-slen # >=0
  if(left)
    return(sprintf("%s%s",substring(padstring,1,padlen),text))
  else
    return(sprintf("%s%s)",text,substring(padstring,1,padlen)))
}


# ExpressionSet: featureData
readExpSet<-function(bam,idx,val="rpmg",phenoData,expData)
{
  # bam: vector of bam-file names
  # idx: optional vector of index-file names
  # val: gptm or rpmg
  # probes: vector containing probe identifier
  # require("Biobase")
  vals<-c("gptm","rpmg")
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
  # Check incoming args
  
  # + + + + + + + + #
  if(!is.character(bam))
    stop("'bam' argument must be character!")
  
  # + + + + + + + + #
  if(!missing(idx))
  {
    if(length(idx)!=length(bam))
      stop("idx must have same length as bam!")
  }
  
  # + + + + + + + + #
  if(!is(phenoData,"AnnotatedDataFrame"))
    stop("[readExpSet] phenoData must be of class 'AnnotatedDataFrame")
  if(nrow(pData(phenoData))!=length(bam))
    stop("[readExpSet]Number of rows in pData(phenoData) must equal length(bam)!")
 
  # + + + + + + + + #
  if(!is.character(val))
    stop("[readExpSet] 'val' argument must be character!")
  if(length(val)>1)
    stop("[readExpSet] 'val' argument must have length 1!")
  if(!is.element(val,vals))
    stop("[readExpSet] 'val' argument must be 'gptm' or 'rpmg'!")
  if(val=="rpmg")
    gptm<-FALSE
  else
    gptm<-TRUE
  
  # + + + + + + + + #
  if(!missing(expData))
  {
    if(!is(expData,"MIAME"))
      stop("[readExpSet] expData must be MIAME")
  }
  
  bm<-Sys.localeconv()[7]
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
  # Read count data from bam-files
  n<-length(bam)
  probes<-rownames(pData(phenoData))

  # + + + + + + + + + + + + + + + + + + #
  # First bam-file
  message("[",format(1,width=2),"/",n,"] ",appendLF=FALSE)
  
  if(missing(idx))
    reader<-bamReader(bam[1],idx=TRUE)
  else
    reader<-bamReader(bam[1],indexname=idx[1])
  
  bsl<-bamGapList(reader)
  dfr<-as.data.frame(bsl)
  bamClose(reader)
  message("size: ",format(size(bsl),big.mark=bm))
  
  if(gptm)
  {
    n_aligns<-nAligns(bsl)
    gptm_fac<-1e7/n_aligns
    dfr$gptm<-round(dfr$nAligns*gptm_fac,5)    
  }else{
    n_gap_aligns<-nAlignGaps(bsl)
    rpmg_fac<-1e6/n_gap_aligns
    dfr$rpmg<-round(dfr$nAligns*rpmg_fac,5)
  }

  res<-dfr[,c("seqid","lstart","lend","rstart","rend",val)]
  # Name new column with probe-id
  names(res)[ncol(res)]<-probes[1]
  # + + + + + + + + + + + + + + + + + + #
  
  for(i in 2:n)
  {
    # + + + + + + + + + + + + + + + + + + #
    # Subsequent bam-file
    message("[",format(i,width=2),"/",n,"] ",appendLF=FALSE)
    
    if(missing(idx))
      reader<-bamReader(bam[i],idx=TRUE)
    else
      reader<-bamReader(bam[i],indexname=idx[i])
    bsl<-bamGapList(reader)
    dfr<-as.data.frame(bsl)
    bamClose(reader)
    message("size: ",format(size(bsl),big.mark=bm))
    
    if(gptm)
    {
      n_aligns<-nAligns(bsl)
      gptm_fac<-1e7/n_aligns
      dfr$gptm<-round(dfr$nAligns*gptm_fac,5)    
    }else{
      n_gap_aligns<-nAlignGaps(bsl)
      rpmg_fac<-1e6/n_gap_aligns
      dfr$rpmg<-round(dfr$nAligns*rpmg_fac,5)
    }  
    
    dfr<-dfr[,c("seqid","lstart","lend","rstart","rend",val)]
    # Name new column with probe-id
    names(dfr)[ncol(dfr)]<-probes[i]
    
    # Merge 
    res<-merge(res,dfr,by=c("seqid","lend","rstart"),all=T,suffixes=c("","2"))
    res$lstart<-pmin(res$lstart,res$lstart2,na.rm=T)
    res$lstart2<-NULL
    res$rend<-pmax(res$rend,res$rend2,na.rm=T)
    res$rend2<-NULL
    # + + + + + + + + + + + + + + + + + + #
  }
  res[is.na(res)]<-0
  res<-res[order(res$seqid,res$lend,res$rstart),]
  
  # Remove possibly prefixing 'chr' (because uninformative)
  seqid<-padd(sub("^chr","",res$seqid),2L,char="0")
  # Unique alphanumerical index
  index<-alphabetIndex(nrow(res),index_alphabet)
  # Index as suffix: fast search (and sort) from left.
  rownames(res)<-paste(index,seqid,sep="")
    
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
  # Create ExpressionSet from data
  
  # assay data
  message("[readExpSet] Create ExpressionSet.")
  ncols<-ncol(res)
  exprm<-as.matrix(res[,6:ncols])
  
  # featureData
  meta<-data.frame(labelDescription=c("sequence identifier","left start","left end", "right start", "right end" ),
                   row.names=c("seqid","lstart","lend","rstart","rend"))
  featureData<-new("AnnotatedDataFrame",data=res[,1:5],varMetadata=meta)
  
  if(missing(expData))
    exSet<-ExpressionSet(assayData=exprm,phenoData=phenoData,featureData=featureData)
  else
    exSet<-ExpressionSet(assayData=exprm,phenoData=phenoData,experimentData=expData,featureData=featureData)
  
  message("[readExpSet] Finished.")
  return(exSet)
}

setGeneric("uniqueJuncAnn",function(object,junc,ann=TRUE,...)standardGeneric("uniqueJuncAnn"))
setMethod("uniqueJuncAnn","ExpressionSet",function(object,junc,ann=TRUE,...){
  
  if(!is(junc,"refJunctions"))
    stop("[uniqueJuncAnn.ExpressionSet] junc must be 'refJunctions'!")
  if(!is.logical(ann))
    stop("[uniqueJuncAnn.ExpressionSet] ann must be logical!")
  
  # Extract sub-tables from ExpressionSet
  message("[uniqueJuncAnn.ExpressionSet] Unifying juncs.")
  uj<-unifyJuncs(junc)
  ex<-exprs(object)
  fd<-featureData(object)@data
  fd$fid<-row.names(fd)
  
  # Merge featureData with annotation data
  message("[uniqueJuncAnn.ExpressionSet] Merge with annotation.")
  mrg<-merge(fd,uj[,c("seqid","lend","rstart","id","gene_id","strand","fexid")],by=c("seqid","lend","rstart"),all.x=TRUE)
  mrg<-mrg[order(mrg$fid),]
  
  # Remove unannotated junction sites
  idna<-is.na(mrg$id)
  bm<-Sys.localeconv()[7]
  if(ann)
  {
    message("[uniqueJuncAnn.ExpressionSet] Remove ",format(sum(idna),big.mark=bm)," unannotated sites.")
    nna<-!idna
    rfd<-mrg[nna,]
    row.names(rfd)<-rfd$fid
    rfd$fid<-NULL
    rex<-ex[nna,]
  }
  else
  {
    message("[uniqueJuncAnn.ExpressionSet] Remove ",format(sum(!idna),big.mark=bm)," annotated sites.")
    rfd<-mrg[idna,]
    row.names(rfd)<-rfd$fid
    rfd$fid<-NULL
    rex<-ex[idna,]
  }
  
  # Assemble result
  message("[uniqueJuncAnn.ExpressionSet] Assemble result.")
  meta<-data.frame(labelDescription=c("sequence identifier","left start","left end", "right start", "right end","id (junction id)","gene id","strand","first exon id" ),
                   row.names=c("seqid","lstart","lend","rstart","rend","id","gene_id","strand","fexid"))
  fdn<-new("AnnotatedDataFrame",data=rfd,varMetadata=meta)
  
  message("[uniqueJuncAnn.ExpressionSet] Finished.")
  return(ExpressionSet(assayData=rex,phenoData=phenoData(object),featureData=fdn))
})

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Read Gene FPKM values from cufflinks files into an ExpressionSet

readCuffGeneFpkm<-function(cuff,phenoData,summ="max")
{
  if(!is.character(cuff))
    stop("[readCuffGeneFpkm] must be character!")
  if(!all(file.exists(cuff)))
    stop("[readCuffGeneFpkm] Non existing files in cuff vector!")
  if(!is(phenoData,"AnnotatedDataFrame"))
    stop("[readCuffGeneFpkm] phenoData must be AnnotatedDataFrame!")
  if(nrow(pData(phenoData))!=length(cuff))
    stop("[readCuffGeneFpkm] Number of rows in pData(phenoData) must equal length(cuff)!")

  if(summ=="max"){
    sm<-1
  }else if(summ=="sum"){
    sm<-2
  }else{
    stop("[readCuffGeneFpkm] summ must be 'max' or 'sum'!")    
  }
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
  # Read count data from bam-files
  n<-length(cuff)
  probes<-rownames(pData(phenoData))
  
  # + + + + + + + + + + + + + + + + + + #
  # First cuff-file
  message("[",format(1,width=2),"/",n,"]\tprobe: ",probes[1])
  tbl<-read.table(cuff[1],sep="\t",header=TRUE)
  #colnames<-c("tracking_id","gene_short_name","locus","FPKM")
  
  # A handful of tracking id's occur up to four times
  # (because of different transcripts)
  # The one with the maximal FPKM is chosen
  track<-unique(tbl$tracking_id)
  mtc<-match(track,tbl$tracking_id)
  gene<-tbl$gene_short_name[mtc]
  locus<-tbl$locus[mtc]
  res<-data.frame(tracking_id=track,gene=gene,locus=locus)
  
  if(sm==1){
    fpkm<-tapply(tbl$FPKM,tbl$tracking_id,max)
  }else if(sm==2){
    fpkm<-tapply(tbl$FPKM,tbl$tracking_id,sum)
  }
  mtc<-match(res$tracking_id,names(fpkm))
  res[,probes[1]]<-fpkm[mtc]
  
  for(i in 2:n)
  {
    # + + + + + + + + + + + + + + + + + + #
    # Subsequent cuff-file
    message("[",format(i,width=2),"/",n,"]\tprobe: ",probes[i])
    tbl<-read.table(cuff[i],sep="\t",header=TRUE)
    
    if(sm==1){
      fpkm<-tapply(tbl$FPKM,tbl$tracking_id,max)
    }else if(sm==2){
      fpkm<-tapply(tbl$FPKM,tbl$tracking_id,sum)
    }  
    
    mtc<-match(names(fpkm),res$tracking_id)
    res[,probes[i]]<-0L
    res[,probes[i]][mtc]<-fpkm
  }
  
  # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
  # Create ExpressionSet from data
  row.names(res)<-res$tracking_id
  exprm<-as.matrix(res[,-(1:3)])
  row.names(exprm)<-res$tracking_id
  
  meta<-data.frame(labelDescription=c("Ensembl gene id", "Gene name","Gene locus"),row.names=c("tracking_id","gene","locus"))
  fd<-new("AnnotatedDataFrame",data=res[,1:3],varMetadata=meta)
  
  exSet<-ExpressionSet(assayData=exprm,phenoData=phenoData,featureData=fd)
  return(invisible(exSet))
}

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  SpliceCountSet                                                                                   #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

#setClass("SpliceCountSet",contains="ExpressionSet")

setMethod("initialize","SpliceCountSet",function(.Object,assayData,phenoData,featureData,exprs,...){
  .Object<-callNextMethod(.Object,phenoData=phenoData,featureData=featureData,exprs=exprs,...) 
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  maxEnt                                                                                           #
#  Calculates splice site scores                                                                    #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

load.maxEnt<-function(file)           
{
  if(missing(file))
    file<-system.file("extdata","maxent.RData",package="spliceSites")
  if(dirname(file)==".")
    mxe<-new("maxEnt",ev=new.env(),basedir=getwd())
  else
    mxe<-new("maxEnt",ev=new.env(),basedir=dirname(file))
  
  load(file,envir=mxe@ev)
  return(invisible(mxe))
}           

read.maxEnt<-function(basedir=getwd())
{
  me<-new("maxEnt",ev=new.env(),basedir=basedir)
  me@ev$me2x5<-scan(file.path(basedir,"me2x5"))
  # scan max-ent files
  maxent_files=file.path(basedir,"splicemodels",c('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9'))
  wmm_files=file.path(basedir,"splicemodels",c('me1s0acc1','me1s0acc2','me1s0acc3','me1s0acc4','me1s0acc5','me1s0acc6','me1s0acc7','me1s0acc8','me1s0acc9'))
  emm_files=file.path(basedir,"splicemodels",c('me2s0acc1','me2s0acc2','me2s0acc3','me2s0acc4','me2s0acc5','me2s0acc6','me2s0acc7','me2s0acc8','me2s0acc9'))  
  for(i in 1:length(maxent_files))
    me@ev$maxent[[i]]<-scan(maxent_files[i])
  for(i in 1:length(wmm_files))
    me@ev$wmm[[i]]<-scan(wmm_files[i])
  for(i in 1:length(emm_files))
    me@ev$emm[[i]]<-scan(emm_files[i])
  return(me)
}

#setGeneric("saveMaxEnt",function(object,filename,...)standardGeneric("saveMaxEnt"))
setMethod("saveMaxEnt","maxEnt",function(object,filename,useBasedir=FALSE,...)
{
  if(!is.character(filename))
    stop("[savemaxEnt.maxEnt] filename must be character!")
  
  if(useBasedir && length(object@basedir)>0)
    save(file=file.path(object@basedir,filename),list=ls(envir=object@ev),envir=object@ev,...)
  else
    save(file=filename,list=ls(envir=object@ev),envir=object@ev,...)
  return(invisible())
})

# Generics from refGenome
#setGeneric("basedir",function(object)standardGeneric("basedir"))
setMethod("basedir","maxEnt",function(object) {return(object@basedir)})
#setGeneric("basedir<-",function(object,value)standardGeneric("basedir<-"))
setReplaceMethod("basedir","maxEnt",function(object,value)
{
  if(!file.exists(value))
    warning("[basedir.maxEnt] Directory '",value,"' does not exist!\n",sep="")
  object@basedir<-value
  return(object)
})

# setGeneric("score5",function(x,seq,pos,...)standardGeneric("score5"))
setMethod("score5","maxEnt",function(x,seq,pos,...){
  if(!is.numeric(pos))
    stop("[score5.maxEnt] pos must be numeric!")
  if(is(seq,"DNAStringSet"))
    seq<-as.character(seq)
  if(!is(seq,"character"))
    stop("[score5.maxEnt] seq must be character (or DNAStringSet)!")

  pos<-as.integer(pos)
  if(length(pos)==1)
    pos<-rep(pos,length(seq))
  if(length(seq)!=length(pos))
    stop("[score5.maxEnt] seq and pos must have same length!")  
  if(any(pos<3))
    stop("[score5.maxEnt] pos must be >=3! At least 3 exon nucs needed!\n")
  if(any(nchar(seq)<pos+6))
    stop("[score5.maxEnt] Sequence length must be >= pos+6=",pos+6,". At least 6 intron nucs needed.\n")
  
  return(.Call("maxent_score5",seq,pos,x@ev$me2x5,PACKAGE="spliceSites"))
})

#setGeneric("scoreSeq5",function(x,seq,frame)standardGeneric("scoreSeq5"))
setMethod("scoreSeq5","maxEnt",function(x,seq,frame){
  if(!is.character(seq))
    stop("[scoreSeq5.maxEnt] seq must be character!")
  
  nc<-nchar(seq)
  if(nc<9) # =l5seq
    stop("[scoreSeq5.maxEnt] nchar(seq) must be >=9!")
  if(missing(frame))
    frame<-c(3,nc-6)
  if(!is.numeric(frame))
    stop("[scoreSeq5.maxEnt] frame must be numeric!")
  if(length(frame)!=2)
    stop("[scoreSeq5.maxEnt] frame must have length 2!")

  #    pos=1-based last exon nuc = 3, l5seq=9
  #    |
  #  ATG | GTC | ATCGAA
  #  123   456   789012
  #    |     |
  #    |      end=6
  #    start=3
  
  return(.Call("maxent_seq_score5",seq,as.integer(frame),x@ev$me2x5,PACKAGE="spliceSites"))
})

setMethod("scoreSeq3","maxEnt",function(x,seq,frame,which="ent",...){
  if(!is.character(seq))
    stop("[scoreSeq3.maxEnt] seq must be character!")
  
  nc<-nchar(seq)
  if(nc<23) # =l3seq
    stop("[scoreSeq3.maxEnt] nchar(seq) must be >=23!")
  if(missing(frame))
    frame<-c(20,nc-3)
  if(!is.numeric(frame))
    stop("[scoreSeq3.maxEnt] frame must be numeric!")
  if(length(frame)!=2)
    stop("[scoreSeq3.maxEnt] frame must have length 2!")
  
  if(which=="ent")
    return(.Call("maxent_seq_score3",seq,as.integer(frame),x@ev$maxent,PACKAGE="spliceSites"))
  if(which=="wmm")
    return(.Call("maxent_seq_score3",seq,as.integer(frame),x@ev$wmm,PACKAGE="spliceSites"))
  if(which=="emm")
    return(.Call("maxent_seq_score3",seq,as.integer(frame),x@ev$emm,PACKAGE="spliceSites"))
})

setMethod("score3","maxEnt",function(x,seq,pos,which="ent",...){
  if(!is.numeric(pos))
    stop("[score3.maxEnt] pos must be numeric!")
  
  if(is(seq,"DNAStringSet"))
    seq<-as.character(seq)
  if(!is(seq,"character"))
    stop("[score3.maxEnt] seq must be character (or DNAStringSet)!")
  
  pos<-as.integer(pos)
  if(length(pos)==1)
    pos<-rep(pos,length(seq))
  if(length(pos)!=length(seq))
    stop("[score3.maxEnt] seq and pos must have same length!")
  if(any(pos<20))
    stop("[score3.maxEnt] pos must be >=20: At least 20 intron nucs needed!\n")
  if(any(nchar(seq)<pos+3))
    stop("[score3.maxEnt] Sequence length must be >= pos+3=",pos+3,". At least 3 exon nucs needed!\n")
  if(!is.character(which))
    stop("[score3.maxEnt] Which argument must be character!")
  if(which=="ent")
    return(.Call("maxent_score3",seq,pos,x@ev$maxent,PACKAGE="spliceSites"))
  if(which=="wmm")
    return(.Call("maxent_score3",seq,pos,x@ev$wmm,PACKAGE="spliceSites"))
  if(which=="emm")
    return(.Call("maxent_score3",seq,pos,x@ev$emm,PACKAGE="spliceSites"))
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  Calculate maxEnts for gapSites objects

#setGeneric("addMaxEnt",function(x,dna,maxent,digits=1)standardGeneric("addMaxEnt"))
setMethod("addMaxEnt","gapSites",function(x,dna,maxent,digits=1){
  if(!is(dna,"DNAStringSet"))
    stop("[addMaxEnd.gapSites] dna must be DNAStringSet!")
  if(!is(maxent,"maxEnt"))
    stop("[addMaxEnd.gapSites] maxEnt must be of class 'maxEnt'")
  if(!is.numeric(digits))
    stop("[addMaxEnd.gapSites] digits must be numeric!")
  if(length(digits)>1)
    stop("[addMaxEnd.gapSites] digits mus have length 1!")
  
  featlen<-3 # nr of exonic   nucleotides
  gaplen<-20 # nr of intronic nucleotides  
  res<-new("gapSites")
  res@dt<-x@dt
  res@nAligns<-x@nAligns
  res@nAlignGaps<-x@nAlignGaps
  
  # Calculate scores for "+"-strand
  lj<-lJunc(x,featlen,gaplen,strand="+")
  dlj<-dnaRanges(lj,dna,useStrand=TRUE,verbose=FALSE)
  
  mtc<-match(res@dt$id,dlj@dt$id)
  if(any(is.na(mtc)))
  {
    bm<-Sys.localeconv()[7]
    message("[addMaxEnt.gapSites] Removing ",format(sum(is.na(mtc)),big.mark=bm),
                      " records due to missing match in dna (DNAStringSet)." )
    res@dt<-x@dt[!is.na(mtc),]
    res@dt$seqid<-factor(res@dt$seqid)
  }
    
  res@dt$mxe_ps5<-round(score5(maxent,dlj@seq,featlen),digits)
  
  rj<-rJunc(x,featlen,gaplen,strand="+",unique=FALSE)
  drj<-dnaRanges(rj,dna,useStrand=TRUE,verbose=FALSE)
  res@dt$mxe_ps3<-round(score3(maxent,drj@seq,gaplen),digits)
  
  # Calculate scores for "-"-strand
  rj<-rJunc(x,featlen,gaplen,strand="-",unique=FALSE)
  drj<-dnaRanges(rj,dna,useStrand=TRUE,verbose=FALSE)
  res@dt$mxe_ms5<-round(score5(maxent,drj@seq,featlen),digits)
  
  lj<-lJunc(x,featlen,gaplen,strand="-",unique=FALSE)
  dlj<-dnaRanges(lj,dna,useStrand=TRUE,verbose=FALSE)
  res@dt$mxe_ms3<-round(score3(maxent,dlj@seq,gaplen),digits)
  
  # Combine score information to strand information
  res@dt$s5strand<-.Call("maxent_score2strand",res@dt$mxe_ps3,res@dt$mxe_ms3)
  res@dt$s3strand<-.Call("maxent_score2strand",res@dt$mxe_ps3,res@dt$mxe_ms3)
  res@dt$meStrand<-.Call("maxent_combine_strand",res@dt$s5strand,res@dt$s3strand)
    
  return(res)
})


#setGeneric("setMeStrand",function(x)standardGeneric("setMeStrand"))
setMethod("setMeStrand","gapSites",function(x){
  mtc<-match("meStrand",names(x@dt))
  if(is.na(mtc))
    stop("[setMeStrand.gapSites] meStrand column not found. Use 'addMaxEnt!")
  
  res<-new("gapSites")
  res@dt<-x@dt
  res@dt$strand<-res@dt$meStrand
  return(res)
})

setMethod("getMeStrand","gapSites",function(x){
  mtc<-match("meStrand",names(x@dt))
  if(is.na(mtc))
    stop("[getMeStrand.gapSites] meStrand column not found. Use 'addMaxEnt!")
  return(x@dt$meStrand)
})

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#  hbond                                                                                            #
#  Calculates splice site scores                                                                    #
#                                                                                                   #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

load.hbond<-function(file)
{
  if(missing(file))
    file<-system.file("extdata","hbs.RData",package="spliceSites")
  
  if(dirname(file)==".")
    hb<-new("hbond",ev=new.env(),basedir=getwd())
  else
    hb<-new("hbond",ev=new.env(),basedir=dirname(file))
  load(file,envir=hb@ev)
  return(invisible(hb))
}

# Generics from refGenome
setMethod("basedir","hbond",function(object) {return(object@basedir)})
setReplaceMethod("basedir","hbond",function(object,value)
{
  if(!file.exists(value))
    warning("[basedir.hbond] Directory '",value,"' does not exist!\n",sep="")
  object@basedir<-value
  return(object)
})

setMethod("hbond","hbond",function(x,seq,pos,...){
  if(!is.numeric(pos))
    stop("[hbond.hbond] pos must be numeric!")
  
  if(is(seq,"DNAStringSet"))
    seq<-as.character(seq)
  if(!is(seq,"character"))
    stop("[hbond.hbond] seq must be character (or DNAStringSet)!")
  
  pos<-as.integer(pos)
  if(length(pos)==1)
    pos<-rep(pos,length(seq))
  if(length(pos)!=length(seq))
    stop("[hbond.hbond] seq and pos must have same length!")
  if(any(pos<3))
    stop("[hbond.hbond] pos must be >=3: 3 exon nucs required!\n")
  if(any(nchar(seq)<pos+8))
    stop("[hbond.hbond] Sequence length must be >= pos+8. At least 8 intron nucs required!\n")
 
  return(.Call("hbond_score",substr(seq,pos-2,pos+8),x@ev$hbs,PACKAGE="spliceSites"))
})

        
#setGeneric("addHbond",function(x,dna,hbond)standardGeneric("addHbond"))
setMethod("addHbond","cdRanges",function(x){
  hb<-load.hbond()
  x@dt$hbond<-hbond(hb,x@seq,x@dt$pos)
  return(x)
})


setMethod("addHbond","gapSites",function(x,dna){
  if(!is(dna,"DNAStringSet"))
    stop("[addMaxEnd.gapSites] dna must be DNAStringSet!")
  
  featlen<-3 # nr of exonic   nucleotides
  gaplen<-8  # nr of intronic nucleotides  
  res<-new("gapSites")
  res@dt<-x@dt
  res@nAligns<-x@nAligns
  res@nAlignGaps<-x@nAlignGaps
  
  # Calculate scores for "+"-strand
  lj<-lJunc(x,featlen,gaplen,strand="+")
  dlj<-dnaRanges(lj,dna,useStrand=TRUE,verbose=FALSE)
  
  mtc<-match(res@dt$id,dlj@dt$id)
  if(any(is.na(mtc)))
  {
    bm<-Sys.localeconv()[7]
    message("[addMaxEnt.gapSites] Removing ",format(sum(is.na(mtc)),big.mark=bm),
                    " records due to missing match in dna (DNAStringSet)." )
    res@dt<-x@dt[!is.na(mtc),]
    res@dt$seqid<-factor(res@dt$seqid)
  } 
  # Loading takes 8 ms so it's ok not to reuse object
  hb<-load.hbond()
  res@dt$lhbond<-.Call("hbond_score",as.character(dlj@seq),hb@ev$hbs,PACKAGE="spliceSites")
  
  rj<-rJunc(x,featlen,gaplen,strand="-")
  drj<-dnaRanges(rj,dna,useStrand=TRUE,verbose=FALSE)
  
  # There should be no mismatches possible
  res@dt$rhbond<-.Call("hbond_score",as.character(drj@seq),hb@ev$hbs,PACKAGE="spliceSites")
  return(res)
})


get_dna_nmers<-function(len)
{
  if(!is.numeric(len))
    stop("[get_dna_nmers] len must be numeric!")
  len<-as.integer(len)
  if(len<=0)
    stop("[get_dna_nmers] len must be >0")
  if(len>32)
    stop("[get_dna_nmers] len must be <=32")
  return(.Call("create_dna_n_mers",len,PACKAGE="spliceSites"))
}


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Plot align depth for splice-junctions
juncplot<-function()
{
  n<-10
  plot(1:n,type="n",axes=FALSE,ann=FALSE)
  op<-par(family="Courier New", lheight=1.5)
  text(1:n,1,LETTERS[1:n])
  par(op)
}

