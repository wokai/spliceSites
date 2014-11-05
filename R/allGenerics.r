
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  File   : allClasses.r                                                        #
#  Date   : 28.Okt.2012                                                         #
#  Content: Generics declarations for package spliceSites                       #
#  Version: 0.7.0                                                               #
#  Author : W. Kaisers                                                          #
#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setGeneric("seqid",           function(x)standardGeneric("seqid"))
setGeneric("start",           function(x)standardGeneric("start"))
setGeneric("end",             function(x)standardGeneric("end"))
setGeneric("id",              function(x)standardGeneric("id"))
setGeneric("count",           function(x)standardGeneric("count"))
setGeneric("seqnames",        function(x)standardGeneric("seqnames"))

setGeneric("lCodons",         function(x,frame=1,keepStrand=TRUE)standardGeneric("lCodons"))
setGeneric("rCodons",         function(x,frame=1L,keepStrand=TRUE)standardGeneric("rCodons"))
setGeneric("lrCodons",        function(x,frame=1L,strand="+")standardGeneric("lrCodons"))
setGeneric("dnaRanges",       function(x,dnaset,useStrand=TRUE,removeUnknownStrand=TRUE,verbose=TRUE,...)standardGeneric("dnaRanges"))

setGeneric("truncateSeq",     function(x,rme=TRUE,trunc=42L)standardGeneric("truncateSeq"))
setGeneric("trypsinCleave",   function(x,minLen=5,...)standardGeneric("trypsinCleave"))


setGeneric("dnaGapSites",     function(x,dnaset,featlen=3,gaplen=8,strand)standardGeneric("dnaGapSites"))
setGeneric("extractRange",    function(object,seqid,start,end)standardGeneric("extractRange"))
setGeneric("seqlogo",         function(x,strand="+",useStrand=FALSE,...) standardGeneric("seqlogo"))

# Basic accessors
setGeneric("seqid" ,          function(x)standardGeneric("seqid"))
setGeneric("lstart",          function(x)standardGeneric("lstart"))
setGeneric("lend"  ,          function(x)standardGeneric("lend"))
setGeneric("rstart",          function(x)standardGeneric("rstart"))
setGeneric("rend"  ,          function(x)standardGeneric("rend"))
setGeneric("gptm",            function(x)standardGeneric("gptm"))
setGeneric("rpmg",            function(x)standardGeneric("rpmg"))
setGeneric("getProfile",      function(x)standardGeneric("getProfile"))
setGeneric("getSequence",     function(x)standardGeneric("getSequence"))

setGeneric("trim_left" ,      function(x,maxlen)standardGeneric("trim_left"))
setGeneric("trim_right",      function(x,maxlen)standardGeneric("trim_right"))
setGeneric("resize_left",     function(x,len)standardGeneric("resize_left"))
setGeneric("resize_right",    function(x,len)standardGeneric("resize_right"))

# Alternative sites
setGeneric("alt_left_ranks" ,  function(x)standardGeneric("alt_left_ranks"))
setGeneric("alt_right_ranks" , function(x)standardGeneric("alt_right_ranks"))

setGeneric("alt_ranks",       function(object)standardGeneric("alt_ranks"))

# Junction Sites
setGeneric("lJunc",          function(x,featlen,gaplen,unique=FALSE,strand,...)standardGeneric("lJunc"))
setGeneric("rJunc",          function(x,featlen,gaplen,unique=FALSE,strand,...)standardGeneric("rJunc"))
setGeneric("lrJunc",         function(x,lfeatlen,rfeatlen,strand,...)standardGeneric("lrJunc"))
setGeneric("lJuncStrand",    function(x,featlen,gaplen,...)standardGeneric("lJuncStrand"))
setGeneric("rJuncStrand",    function(x,featlen,gaplen,...)standardGeneric("rJuncStrand"))
setGeneric("lrJuncStrand",   function(x,lfeatlen,rfeatlen,...)standardGeneric("lrJuncStrand"))

setGeneric("annotate",       function(object,genome)standardGeneric("annotate"))
setGeneric("getLann" ,       function(x)standardGeneric("getLann"))
setGeneric("getRann" ,       function(x)standardGeneric("getRann"))
setGeneric("getAnnStrand",   function(object)standardGeneric("getAnnStrand"))
setGeneric("addGeneAlignPart",function(x)standardGeneric("addGeneAlignPart"))
setGeneric("sortTable",      function(x)standardGeneric("sortTable"))


# write.files
setGeneric("write.files",    function(x,path,filename,...)standardGeneric("write.files"))
setGeneric("write.annDNA.tables",function(object,dnaset,filename,featlen=3,gaplen=8,sep=";",dec=".",row.names=FALSE)standardGeneric("write.annDNA.tables"))


# plot Distance between lend and left_end and between rstart and right_start
setGeneric("plot_diff",		function(x,n=20)standardGeneric("plot_diff"))
setGeneric("plot_diff_ranks",	function(x)standardGeneric("plot_diff_ranks"))

# Generics for keyProfiler
setGeneric("addKeyTable",	function(object,keyTable,index,values)standardGeneric("addKeyTable"))
setGeneric("getKeyTable",	function(object)standardGeneric("getKeyTable"))
setGeneric("appendKeyTable",	function(object,keytable,prefix,valFactor,rateFactor,digits)standardGeneric("appendKeyTable"))

# Generics for maxEntScore
setGeneric("saveMaxEnt",	function(object,filename,...)standardGeneric("saveMaxEnt"))
setGeneric("score5",		function(x,seq,pos,...)standardGeneric("score5"))
setGeneric("score3",		function(x,seq,pos,which="ent",...)standardGeneric("score3"))
setGeneric("addMaxEnt",		function(x,dna,maxent,digits=1)standardGeneric("addMaxEnt"))
setGeneric("setMeStrand",	function(x)standardGeneric("setMeStrand"))
setGeneric("getMeStrand",	function(x)standardGeneric("getMeStrand"))
setGeneric("scoreSeq5",		function(x,seq,frame)standardGeneric("scoreSeq5"))
setGeneric("scoreSeq3",		function(x,seq,frame,which="ent",...)standardGeneric("scoreSeq3"))

# Generics for hbond Score
setGeneric("hbond",		function(x,seq,pos,...)standardGeneric("hbond"))
setGeneric("addHbond",		function(x,dna)standardGeneric("addHbond"))



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Generics provided by packages on which spliceSites depends on:

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Provided by IRanges >=1.18.2
# setGeneric("width",           function(x)standardGeneric("width"))

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Provided by Biostrings >=2.28.0 
# setGeneric("translate",      function(x)standardGeneric("translate"))

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Provided by BiocGenerics >=0.6.0
# setGeneric("annotation",   function(object)standardGeneric("annotation"))
# setGeneric("annotation<-", function(object,value)standardGeneric("annotation<-"))
# setGeneric("strand",       function(x,...)standardGeneric("strand"))

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Provided by rbamtools (specialized for gapSites)
#setGeneric("nAligns",         function(object)standardGeneric("nAligns"))
#setGeneric("nAlignGaps",      function(object)standardGeneric("nAlignGaps"))

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Provided by refGenome (specialized for maxEntScore)

#setGeneric("basedir",function(object)standardGeneric("basedir"))
#setGeneric("basedir<-",function(object,value)standardGeneric("basedir<-"))
#setGeneric("extractByGeneName",function(object,geneNames,src,...)standardGeneric("extractByGeneName"))


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
