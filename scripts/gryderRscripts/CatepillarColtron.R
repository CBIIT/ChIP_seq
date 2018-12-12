
setwd("K:/projects/ChIP_seq")
project.folder = "projects/MCC/coltron/"

###PRELOAD dataframe "coltron" and character value "SampleList" from program gatherColtron.R

#make coltron derived dataframe with $from and $to with $width
#get expression data to rank TFs 
#make coltron derived $vsize based on expression
#make Network file from Edge list 

library(qgraph)
project.prefix = "biggernodecoral2"
edgecolor = "coral2"

lapply(SampleList, function(x) {
  filepath <- subset(coltron, coltron$SampleFiles %in% x) #isolate single sample
  library(plyr)
  
  #load SE TFs and rank by expression
    TF_Enhancers = read.table(paste(filepath$path,filepath$SampleFiles,"_ENHANCER_ASSIGNMENT.txt",sep=""),sep="\t", header=F)
    TF_Expression = TF_Enhancers[,c(1,2,3)]
    colnames(TF_Expression) = c("RefseqID","GeneID","FPKM")
    TF_Expression = TF_Expression[!duplicated(TF_Expression[,c('GeneID')]),]  #remove if GeneID has been duplicated
    TF_Expression = TF_Expression[order(-TF_Expression$FPKM), ]  #sort by FPKM
    TF_Expression$rank.order<-as.numeric(factor(TF_Expression$FPKM,levels=TF_Expression$FPKM)) #give rank order
  #load TF edges and add rank to edges
  TF_Edges <- read.table(paste(filepath$path,filepath$SampleFiles,"_EDGE_LIST.txt",sep=""),sep="\t", header=T)
    colnames(TF_Edges) = c("GeneID","To")
    TF_Edges <- join(TF_Edges, TF_Expression, by = "GeneID")
    colnames(TF_Edges) = c("From","GeneID","FromRef","FromFPKM","FromRank")
    TF_Edges <- join(TF_Edges, TF_Expression, by = "GeneID")
    colnames(TF_Edges) = c("From","To","FromRef","FromFPKM","FromRank","ToRef","ToFPKM","ToRank")
    TF_Edges$Fromlog2FPKM = log(TF_Edges$FromFPKM,2)
      write.table(TF_Edges,paste(project.folder,project.prefix,x,"TF_Edges.txt",sep=""),sep="\t",row.names=F, quote=F, col.names=T)
    TF_CatepillarEdges = TF_Edges[,c("FromRank","ToRank","Fromlog2FPKM")] #Catepillar dataframe, with edge weight based on expression of "from" column
  #make vsize dataframe based on nodelist
  TF_Nodes <- as.data.frame(read.table(paste(filepath$path,filepath$SampleFiles,"_NODELIST.txt",sep=""),sep="\t", header=F))
    colnames(TF_Nodes) = c("GeneID")
  TF_Nodes.ranked = join(TF_Nodes, TF_Expression, by = "GeneID")
    TF_CatepillarSize = TF_Nodes.ranked[,c("rank.order","FPKM")]
    TF_CatepillarSize$log2FPKM = log(TF_CatepillarSize$FPKM,2)
    TF_CatepillarSize = TF_CatepillarSize[order(-TF_CatepillarSize$FPKM), ]  #sort by FPKM
    CatepillarLayout = matrix(1:nrow(TF_CatepillarSize),nrow=1)
    
    CatepillarLayoutSpaced = TF_CatepillarSize[,c(1,3)]
      CatepillarLayoutSpaced$Xcoord = 0
      CatepillarLayoutSpaced$Ycoord = 12/(CatepillarLayoutSpaced$rank.order+12)  
      CatepillarLayoutSpaced = as.matrix(CatepillarLayoutSpaced[,c(3,4)])
  pdf(paste(project.folder,project.prefix,x,"_Catepillar_All_SE_TFs.pdf",sep=""))
    qgraph(TF_CatepillarEdges, mode = "strength", minimum = 2, esize = 2, edge.color = edgecolor, asize=2, labels=TF_Expression$GeneID, vsize = TF_CatepillarSize$log2FPKM/2, layout=CatepillarLayoutSpaced, layoutScale=c(0.5,1),directed = T)
    mtext(paste(project.prefix,x,"_Catepillar_All_SE_TFs",sep=""), side=1)
    dev.off()
  pdf(paste(project.folder,project.prefix,x,"_Spider_All_SE_TFs.pdf",sep=""))
    qgraph(TF_CatepillarEdges, mode = "strength", minimum = 2, esize = 2, edge.color = edgecolor, asize=2, labels=TF_Expression$GeneID, vsize = TF_CatepillarSize$log2FPKM/1, layoutScale=c(1,0.7),directed = T)
    mtext(paste(project.prefix,x,"_Spider_All_SE_TFs",sep=""), side=1)
    dev.off()
    
  #load degree table and remove TFs which have no Out degree, then recreate plots
  TF_Degree <- read.table(paste(filepath$path,filepath$SampleFiles,"_DEGREE_TABLE.txt",sep=""),sep="\t", header=T)
    TF_Out = subset(TF_Degree, TF_Degree$Out_Degree>0)
      
    #load SE TFs and rank by expression
 
    TF_ExpressionOut = subset(TF_Expression,TF_Expression$GeneID %in% TF_Out$Tf)  #remove TFs lacking OUT binding
      TF_ExpressionOut = TF_ExpressionOut[,c(1,2,3)]
      TF_ExpressionOut = TF_ExpressionOut[order(-TF_ExpressionOut$FPKM), ]  #sort by FPKM
      TF_ExpressionOut$rank.order<-as.numeric(factor(TF_ExpressionOut$FPKM,levels=TF_ExpressionOut$FPKM)) #give rank order
      
    TF_EdgesOut = subset(TF_Edges,TF_Edges$From %in% TF_Out$Tf)
      TF_EdgesOut = subset(TF_EdgesOut,TF_EdgesOut$To %in% TF_Out$Tf)
      #deleate old ranks and add new ranks to edges
        TF_EdgesOut = TF_EdgesOut[,c("From","To")]
        colnames(TF_EdgesOut) = c("GeneID","To")
        TF_EdgesOut <- join(TF_EdgesOut, TF_ExpressionOut, by = "GeneID")
        colnames(TF_EdgesOut) = c("From","GeneID","FromRef","FromFPKM","FromRank")
        TF_EdgesOut <- join(TF_EdgesOut, TF_ExpressionOut, by = "GeneID")
        colnames(TF_EdgesOut) = c("From","To","FromRef","FromFPKM","FromRank","ToRef","ToFPKM","ToRank")
        TF_EdgesOut$Fromlog2FPKM = log(TF_EdgesOut$FromFPKM,2)
      TF_CatepillarEdgesOut = TF_EdgesOut[,c("FromRank","ToRank","Fromlog2FPKM")]
      
      #make vsize dataframe based on nodelist
      TF_NodesOut <- as.data.frame(read.table(paste(filepath$path,filepath$SampleFiles,"_NODELIST.txt",sep=""),sep="\t", header=F))
        colnames(TF_NodesOut) = c("GeneID")
        TF_NodesOut = subset(TF_NodesOut,TF_NodesOut$GeneID %in% TF_Out$Tf)  #remove TFs lacking OUT binding
        TF_NodesOut.ranked = join(TF_NodesOut, TF_ExpressionOut, by = "GeneID")
      
      TF_CatepillarSizeOut = TF_NodesOut.ranked[,c("rank.order","FPKM")]
        TF_CatepillarSizeOut$log2FPKM = log(TF_CatepillarSizeOut$FPKM,2)
        TF_CatepillarSizeOut = TF_CatepillarSizeOut[order(-TF_CatepillarSizeOut$FPKM), ]  #sort by FPKM
      
      CatepillarLayoutOut = matrix(1:nrow(TF_CatepillarSizeOut),nrow=1)
        CatepillarLayoutSpacedOut = TF_CatepillarSizeOut[,c(1,3)]
        CatepillarLayoutSpacedOut$Xcoord = 0
        CatepillarLayoutSpacedOut$Ycoord = 12/(CatepillarLayoutSpacedOut$rank.order+12)
        CatepillarLayoutSpacedOut = as.matrix(CatepillarLayoutSpacedOut[,c(3,4)])
      
      pdf(paste(project.folder,project.prefix,x,"_Catepillar_Out_SE_TFs.pdf",sep=""))
        qgraph(TF_CatepillarEdgesOut, mode = "strength", minimum = 2, esize = 2, edge.color = edgecolor, asize=2, labels=TF_ExpressionOut$GeneID, vsize = TF_CatepillarSizeOut$log2FPKM/2, layout=CatepillarLayoutSpacedOut, layoutScale=c(0.5,1),directed = T)
        mtext(paste(project.prefix,x,"_Catepillar_Out_SE_TFs",sep=""), side=1)
        dev.off()
      pdf(paste(project.folder,project.prefix,x,"_Spider_Out_SE_TFs.pdf",sep=""))
        qgraph(TF_CatepillarEdgesOut, mode = "strength", minimum = 2, esize = 2, edge.color = edgecolor, asize=2, labels=TF_ExpressionOut$GeneID, vsize = TF_CatepillarSizeOut$log2FPKM/1, layoutScale=c(1,0.7),directed = T)
        mtext(paste(project.prefix,x,"_Spider_Out_SE_TFs",sep=""), side=1)
        dev.off()
})


#other random ways of graphing Coltron TF data
  qgraph(TF_CatepillarEdges, mode = "strength", vsize = TF_CatepillarSize$log2FPKM,layout = CatepillarLayout, directed = T)
  library(igraph)
  qgraph(TF_CatepillarEdges, mode = "strength", vsize = TF_CatepillarSize$log2FPKM, layout=layout_with_fr, directed = T)
  qgraph(TF_CatepillarEdges, mode = "strength", vsize = TF_CatepillarSize$log2FPKM, layout="circle" , directed = T)
