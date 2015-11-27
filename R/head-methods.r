
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## head-methods
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("head","gapSites",function(x,n=6L,digits=10,...)
{
    bm<-Sys.localeconv()[7]
    nRows<-nrow(x@dt)
    cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
    
    cat("nAligns: ",format(x@nAligns,big.mark=bm,digits=digits),
        "\tnAlignGaps: ",format(x@nAlignGaps,big.mark=bm,digits=digits),"\n")
    
    if(nRows>0)
    {
        ## Add annotation columns for printout
        res<-x@dt[1:min(n,nRows),]
        if(!is.null(x@annotation))
            res<-merge(res,x@annotation[,c("id","gene_id","gene_name")],by="id")
        print(res)
    }
})


setMethod("head","dnaGapSites",function(x,n=6L,digits=10,...)
{
    bm<-Sys.localeconv()[7]
    nRows<-nrow(x@dt)
    cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
    
    cat("nAligns: ",format(x@nAligns,big.mark=bm,digits=digits),
        "\tnAlignGaps: ",format(x@nAlignGaps,big.mark=bm,digits=digits),"\n")
    
    n<-min(n,nRows)
    
    dt<-x@dt[1:n,]
    dt$seq<-as.character(x@seq[1:n])
    
    if(!is.null(x@annotation))
    {
        dt$gene_id<-x@annotation$gene_id[1:n]
        dt$gene_name<-x@annotation$gene_name[1:n]
    }    
    print(dt)
})


setMethod("head","aaGapSites",function(x,n=6L,digits=10,...)
{
    bm<-Sys.localeconv()[7]
    nRows<-nrow(x@dt)
    cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
    
    cat("nAligns: ",format(x@nAligns,big.mark=bm,digits=digits),
        "\tnAlignGaps: ",format(x@nAlignGaps,big.mark=bm,digits=digits),"\n")
    
    n<-min(n,nRows)
    dt<-x@dt[1:n,]
    dt$seq<-as.character(x@seq[1:n])
    
    if(!is.null(x@annotation))
    {
        dt$gene_id<-x@annotation$gene_id[1:n]
        dt$gene_name<-x@annotation$gene_name[1:n]
    }    
    print(dt)
})

setMethod("head","cRanges",function(x,n=6L,digits=10,...)
{
    bm<-Sys.localeconv()[7]
    nRows<-nrow(x@dt)
    cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
    dt<-x@dt[1:min(nRows,n),]
    
    print(dt)
})

setMethod("head","cdRanges",function(x,n=6L,digits=10,...)
{
    bm<-Sys.localeconv()[7]
    nRows<-nrow(x@dt)
    cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
    dt<-x@dt[1:min(nRows,n),]
    dt$seq<-as.character(x@seq[1:min(nRows,n)])
    
    print(dt)
})

setMethod("head","caRanges",function(x,n=6L,digits=10,...)
{
    bm<-Sys.localeconv()[7]
    nRows<-nrow(x@dt)
    cat("Object of class",class(x),"with",format(nRows,big.mark=bm),"rows and",ncol(x@dt),"columns.\n")
    dt<-x@dt[1:min(nRows,n),]
    dt$seq<-as.character(x@seq[1:min(nRows,n)])
    
    print(dt)
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
