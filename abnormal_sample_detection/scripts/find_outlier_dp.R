# find outliers from linear regression between Mean and SD of region depth
# Input: pooled region depth
# Output: outliers individuals and regression graph

#if (!require(car)) {
    #install.packages("car", repos="http://cran.stat.ucla.edu/")
#}

find.outlier.dp <- function(data, figure.file) {
    require(car)
    fit <- lm(SD ~ Mean, data=data)
    outlier <- names(outlierTest(fit,n.max=10,cutoff=0.05)$p)
    #index <- which(rownames(data) %in% outlier)
    #print(fit)
    png(figure.file)
    plot(data[,c("Mean","SD")])
    points(data[outlier,c("Mean","SD")],col="red")
    abline(fit)
    dev.off()
    #print(cor(data[,c("Mean","SD")]))
    #print(cor(data[-index,c("Mean","SD")]))
    message("\nProblematic samples with abnormal mean and standard deviation of depth: ")
    print(data[outlier,c("Mean","SD")])
}

args <- commandArgs(trailingOnly = TRUE)
file.list <- args[1]
result.file <- args[2]
pngfile <- args[3]
files <- scan(file.list, character())
for (file in files) {
  temp <- read.csv(file, check.names=F)
  if (!exists("region.mean")) {
    region.mean <- temp
  } else {
    region.mean <- rbind(region.mean, temp)
  }
}

region.mean <- region.mean[, 2:455]
result <- data.frame(row.names = colnames(region.mean))
result$Mean <- apply(region.mean, 2, function(col) mean(col[col < 100], na.rm=T))
result$SD <- apply(region.mean, 2, function(col) sd(col, na.rm=T))
result <- result[order(rownames(result)), ]
find.outlier.dp(result, pngfile)
write.table(result, result.file, quote=F)
