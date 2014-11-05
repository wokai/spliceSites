

setMethod("head","gapSites",function(x,n=6L,digits=10,...){
  bm<-Sys.localeconv()[7]
  nRows<-nrow(x@dt)
  cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
  cat("nAligns: ",format(x@nAligns,big.mark=bm,digits=digits),
      "\tnAlignGaps: ",format(x@nAlignGaps,big.mark=bm,digits=digits),"\n")
  
  if(nRows>0)
  {
    # Add annotation columns for printout
    res<-x@dt[1:min(n,nRows),]
    if(!is.null(x@lann))
      res<-merge(res,x@lann[,c("id","left_gene_id","left_gene_name")],by="id")
    if(!is.null(x@rann))
      res<-merge(res,x@rann[,c("id","right_gene_id","right_gene_name")],by="id")
    print(res)    
  }
})

setMethod("head","dnaGapSites",function(x,n=6L,digits=10,...){
  bm<-Sys.localeconv()[7]
  nRows<-nrow(x@dt)
  cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
  cat("nAligns: ",format(x@nAligns,big.mark=bm,digits=digits),
      "\tnAlignGaps: ",format(x@nAlignGaps,big.mark=bm,digits=digits),"\n")
  n<-min(n,nRows)
  dt<-x@dt[1:n,]
  dt$seq<-as.character(x@seq[1:n])
  if(!is.null(x@lann))
  {
    dt$left_gene_id<-x@lann$left_gene_id[1:n]
    dt$left_gene_name<-x@lann$left_gene_name[1:n]
  }
  if(!is.null(x@rann))
  {
    dt$right_gene_id<-x@rann$right_gene_id[1:n]
    dt$right_gene_name<-x@rann$right_gene_name[1:n]
  }
  print(dt)
})

setMethod("head","aaGapSites",function(x,n=6L,digits=10,...){
  bm<-Sys.localeconv()[7]
  nRows<-nrow(x@dt)
  cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
  cat("nAligns: ",format(x@nAligns,big.mark=bm,digits=digits),
      "\tnAlignGaps: ",format(x@nAlignGaps,big.mark=bm,digits=digits),"\n")
  n<-min(n,nRows)
  dt<-x@dt[1:n,]
  dt$seq<-as.character(x@seq[1:n])
  if(!is.null(x@lann))
  {
    dt$left_gene_id<-x@lann$left_gene_id[1:n]
    dt$left_gene_name<-x@lann$left_gene_name[1:n]
  }
  if(!is.null(x@rann))
  {
    dt$right_gene_id<-x@rann$right_gene_id[1:n] 
    dt$right_gene_name<-x@rann$right_gene_name[1:n]
  }
  print(dt)
})

setMethod("head","cRanges",function(x,n=6L,digits=10,...){
  bm<-Sys.localeconv()[7]
  nRows<-nrow(x@dt)
  cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
  dt<-x@dt[1:min(nRows,n),]
  print(dt)
})

setMethod("head","cdRanges",function(x,n=6L,digits=10,...){
  bm<-Sys.localeconv()[7]
  nRows<-nrow(x@dt)
  cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
  dt<-x@dt[1:min(nRows,n),]
  dt$seq<-as.character(x@seq[1:min(nRows,n)])
  print(dt)
})

setMethod("head","caRanges",function(x,n=6L,digits=10,...){
  bm<-Sys.localeconv()[7]
  nRows<-nrow(x@dt)
  cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
  dt<-x@dt[1:min(nRows,n),]
  dt$seq<-as.character(x@seq[1:min(nRows,n)])
  print(dt)
})

