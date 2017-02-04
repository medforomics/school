#!/cm/shared/apps/R/intel/3.2.1/bin/Rscript
library(edgeR)
library(DESeq2)
library("RColorBrewer")
library("gplots")
library(qusage)

rowMax <- function(x) apply(x,1,max)
col.grp <- function(n,b) {
  colorvec <- vector(mode="character", length=length(n))
  plotcols <- rainbow(length(b))
  for (i in 1:length(n)) {
    for (j in 1:length(b)) {
      if ( n[i] == b[j] ) {
        colorvec[i] = plotcols[j]
      }
    }
  }
  c(colorvec)
}


#################### Read in Data ################################
genenames <- read.table(file="genenames.txt",header=TRUE,sep='\t')

tbl <- read.table('countTable.txt',header=TRUE,sep="\t")
tbl2 <- read.table('countTable.logCPM.txt',header=TRUE,sep="\t")
ct <- tbl[,4:length(tbl)]
row.names(ct) <- tbl$ENSEMBL

samples<- names(ct)

dtbl <- read.table('design.txt',header=TRUE,sep="\t")
samtbl <- merge(as.data.frame(samples),dtbl,by.x="samples",by.y="SampleID",all.x=TRUE,sort=FALSE)
grpnames <- levels(factor(as.character(samtbl$SampleGroup)))
samtbl$SampleGroup <- factor(samtbl$SampleGroup, levels=grpnames)

colData <- samtbl[c('SampleGroup','SubjectID')]
row.names(colData) <- samtbl$samples

dds <- DESeqDataSetFromMatrix(countData=ct,colData= colData,design= ~ SampleGroup)
dds <- dds[ rowMax(counts(dds)) > 30, ]
dds <- dds[ colSums(counts(dds)) > 1000000]

countTable <- counts(dds)
grps <- as.character(samtbl$SampleGroup)
col.blocks <-col.grp(grps,levels(factor(grps)))
libSizes <- as.vector(colSums(countTable))

logcpm <- tbl2[,4:length(tbl2)]
row.names(logcpm) <- tbl2$SYMBOL

tmp.tab <- aggregate(t(logcpm),by=list(as.character(samtbl$SampleGroup)),FUN=mean)
row.names(tmp.tab) <- tmp.tab$Group.1
mean.by.group <- round(log2(t(tmp.tab[,c(2:ncol(tmp.tab))])+1), digits = 2)

#################### Run DESEQ2 ################################
dds <- DESeq(dds)
rld <- rlogTransformation(dds, blind=TRUE)
sampleDists <- dist(t(assay(rld)))

png(file="samples_heatmap.png",bg ="transparent",height=768,width=1024)
heatmap.2(as.matrix(sampleDists), col = bluered(100),RowSideColors = col.blocks,srtRow=45,srtCol=45,trace="none", margins=c(5, 5))
dev.off()

#Compare Samples using PCA
png(file="pca.png",bg ="transparent",height=768,width=1024)
print(plotPCA(rld, intgroup="SampleGroup"),col.hab=col.blocks)
dev.off()

#Do all pairwise comparisons
contrast <- resultsNames(dds)
cond<- levels(colData(dds)$SampleGroup)
a <- length(cond)-1
for (i in 1:a) {
  for (j in 2:length(cond)) {
    if (i == j) {
      next
    } else {
      res <- as.data.frame(results(dds,contrast=c("SampleGroup",cond[i],cond[j])))
      res2 <- merge(genenames,res,by.x='ensembl',by.y='row.names',all.y=TRUE,all.x=FALSE)
      output <- merge(res2,mean.by.group,by.y="row.names",by.x='symbol')
      output$rawP <- output$pvalue
      output$logFC <- output$log2FoldChange
      output$fdr <- output$padj
      output$bonf <- p.adjust(output$rawP, method ='bonferroni')
      write.table(output,file=paste(cond[i],'_',cond[j],'.deseq2.txt',sep=""),quote=FALSE,row.names=FALSE,sep='\t')
      filt.out <- na.omit(output[output$fdr < 0.05,])
      if (nrow(filt.out) > 2) {
      	 subset <- logcpm[row.names(logcpm) %in% filt.out$symbol,]
      	 gnames <- filt.out[c('ensembl','symbol')]
      	 s <- merge(gnames,subset,by.x="ensembl",by.y="row.names",all.x=FALSE,all.y=TRUE,sort=FALSE)
      	 STREE <- hclust(dist(t(subset)))
      	 zscores <- scale(t(subset))
      	 ngenes <- length(colnames(zscores))
      	 textscale <- (1/(ngenes/30))
      	 if (textscale > 1) {
            textscale <-1
      	 }
      	 if (textscale < 0.1) {
            textscale <- 0.1
      	 }
      	 png(file=paste(cond[i],'_',cond[j],'.heatmap.deseq2.png',sep=""),height=768,width=1024)
      	 heatmap.2(zscores, col = bluered(100),Rowv = as.dendrogram(STREE), RowSideColors = col.blocks,dendrogram='row', cexCol=textscale,labCol=s$symbol,srtRow=45,srtCol=45,trace="none", margins=c(5, 5))
      	 legend("topright",legend=grpnames,col=rainbow(length(grpnames)),pch=20,cex=0.5)
      	 dev.off()
      }
    }
  }
}

###################### Run EdgeR ########################
###################### Run QuSage ######################

MSIG.geneSets <- read.gmt('geneset.gmt')

design <- model.matrix(~grps)
d <- DGEList(counts=countTable,group=grps,lib.size=libSizes)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
png(file="mds.png",bg ="transparent",height=768,width=1024)
plotMDS(d, labels=grps,col=col.blocks)
legend("topleft",legend=grpnames,col=rainbow(length(grpnames)),pch=20)
dev.off()
cond <-levels(d$samples$group)
colnames(design) <- levels(d$samples$group)
a <- length(cond)-1
for (i in 1:a) {
  for (j in 2:length(cond)) {
    if (i == j) {
      next
    } else {
      c <- exactTest(d, pair=c(cond[j],cond[i]))
      res <- c$table
      res2 <- merge(genenames,res,by.x='ensembl',by.y='row.names',all.y=TRUE,all.x=FALSE)
      output <- merge(res2,mean.by.group,by.y="row.names",by.x='symbol')
      output$rawP <- output$PValue
      output$logFC <- output$logFC
      output$fdr <- p.adjust(output$rawP, method ='fdr')
      output$bonf <- p.adjust(output$rawP, method ='bonferroni')
      write.table(output,file=paste(cond[i],'_',cond[j],'.edgeR.txt',sep=""),quote=FALSE,row.names=FALSE,sep='\t')
      filt.out <- na.omit(output[output$fdr < 0.05,])
      if (nrow(filt.out) > 2) {
      subset <- logcpm[row.names(logcpm) %in% filt.out$symbol,]
      gnames <- filt.out[c('ensembl','symbol')]
      s <- merge(gnames,subset,by.x="ensembl",by.y="row.names",all.x=FALSE,all.y=TRUE,sort=FALSE)
      STREE <- hclust(dist(t(subset)))
      zscores <- scale(t(subset))
      ngenes <- length(colnames(zscores))
      textscale <- (1/(ngenes/30))
      if (textscale > 1) {
        textscale <-1
      }
      if (textscale < 0.1) {
        textscale <- 0.1
      }
      png(file=paste(cond[i],'_',cond[j],'.heatmap.edgeR.png',sep=""),height=768,width=1024)
      heatmap.2(zscores, col = bluered(100),Rowv = as.dendrogram(STREE), RowSideColors = col.blocks,dendrogram='row', cexCol=textscale,labCol=s$symbol,srtRow=45,srtCol=45,trace="none", margins=c(5, 5))
      legend("topright",legend=grpnames,col=rainbow(length(grpnames)),pch=20,cex=0.5)
      dev.off()
      gcont <- paste(cond[j],cond[i],sep='-')
      qs.results = qusage(logcpm, grps,gcont,MSIG.geneSets)
      save(qs.results,file=paste(cond[i],'_',cond[j],'.qusage.rda',sep=""))
      }
      }
  }
}

###################### Run Limma VOOM ########################

# design <- model.matrix(~0+grps)
# colnames(design) <- grpnames
# a <- length(cond)-1
# design.pairs <- c()
# k <- 0
# for (i in 1:a) {
#   for (j in 2:length(cond)) {
#     if (i == j) {
#       next
#     } else {
#     k <- k+1
#     design.pairs[k] <- paste(cond[i],'-',cond[j],sep='')
#     }
#     }
#     }

#contrast.matrix <- makeContrasts(design.pairs,levels=design)

#d <- DGEList(counts=countTable,group=grps,lib.size=libSizes)
#d <- calcNormFactors(d)
#d <- estimateCommonDisp(d)

# v <- voom(d,design,plot=TRUE)
# fit <- lmFit(v,design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)

# comps <- colnames(fit2$coefficients)
# for (i in 1:length(comps)) {
#   res <- topTable(fit2,coef=i,number=Inf,sort.by="P")
#   res2 <- merge(genenames,res,by.x='ensembl',by.y='row.names',all.y=TRUE,all.x=FALSE)
#   output <- merge(res2,mean.by.group,by.y="row.names",by.x='symbol')
#   output$rawP <- output$P.Value
#   output$logFC <- output$adj.P.Val
#   output$fdr <- p.adjust(output$rawP, method ='fdr')
#   output$bonf <- p.adjust(output$rawP, method ='bonferroni')
#   write.table(output,file=paste(cond[i],'_',cond[j],'.voom.txt',sep=""),quote=FALSE,row.names=FALSE,sep='\t')
#   filt.out <- na.omit(output[output$fdr < 0.05,])
#   if (nrow(filt.out) > 2) {
#   subset <- logcpm[row.names(logcpm) %in% filt.out$ensembl,]
#   gnames <- filt.out[c('ensembl','symbol')]
#   s <- merge(gnames,subset,by.x="ensembl",by.y="row.names",all.x=FALSE,all.y=TRUE,sort=FALSE)
#   STREE <- hclust(dist(t(subset)))
#   zscores <- scale(t(subset))
#   ngenes <- length(colnames(zscores))
#   textscale <- (1/(ngenes/30))
#   if (textscale > 1) {
#      textscale <-1
#   }
#   if (textscale < 0.1) {
#     textscale <- 0.1
#   }
#   png(file=paste(cond[i],'_',cond[j],'.heatmap.voom.png',sep=""),height=768,width=1024)
#   heatmap.2(zscores, col = bluered(100),Rowv = as.dendrogram(STREE), RowSideColors = col.blocks,dendrogram='row', cexCol=textscale,labCol=s$symbol,srtRow=45,srtCol=45,trace="none", margins=c(5, 5))
#   legend("topright",legend=grpnames,col=rainbow(length(grpnames)),pch=20,cex=0.5)
#   dev.off()
#   }
# }


