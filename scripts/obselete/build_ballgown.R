library(ballgown)

args<-commandArgs(TRUE)

sample_dirs <- args
samples <- gsub('_stringtie','',sample_dirs)
design <- "design.txt"

samtbl <- read.table(file=design,header=TRUE,sep='\t')

bg <- ballgown(samples=sample_dirs, meas='all')
samples <- gsub('_stringtie','',sampleNames(bg))
mergetbl <- merge(as.data.frame(samples),samtbl,by.x="samples",by.y="SampleID",all.x=TRUE,sort=FALSE)
pData(bg) = data.frame(id=samples, group=as.character(mergetbl$SampleGroup))

save(bg,file="bg.rda")
