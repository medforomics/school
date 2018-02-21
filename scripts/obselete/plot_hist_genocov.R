print(files <- list.files(path='.',pattern="genomecov.txt$"))
print(labs <- gsub(".genomecov.txt", "", files, perl=TRUE))
i <- 1

cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
    tbl <- read.table(files[i])
    bins.x <- cut(tbl$V2,breaks=c(100*-1:20,Inf),labels=c(as.character(100*0:20),'2000+'))
    x <- levels(bins.x)
    y <- aggregate(tbl$V3, by=list(Category=bins.x), FUN=sum)
    png(paste(labs[i],".coverage_histogram.png",sep=''), h=1000, w=1000, pointsize=20)
    barplot(y$x,names.arg=as.character(y$Category),las=2)
    dev.off()
}   

# Pick some colors
# Ugly:
# cols <- 1:length(cov)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
library(RColorBrewer)
cols <- rainbow(length(files))

# Save the graph to a file
png("coverage_cdf.png", h=1000, w=1000, pointsize=20)

# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[[1]][2:1001, 2], cov_cumul[[1]][1:1000], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage")
abline(v = 50, col = "gray60")
abline(v = 100, col = "gray60")
abline(v = 200, col = "gray60")
abline(v = 500, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(50,100,200,500), labels=c(50,100,200,500))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:1001, 2], cov_cumul[[i]][1:1000], type='l', lwd=3, col=cols[i])

# Add a legend using the nice sample labeles rather than the full filenames.
legend("bottom", legend=labs, col=cols, lty=1, lwd=4,cex=0.5)

dev.off()
