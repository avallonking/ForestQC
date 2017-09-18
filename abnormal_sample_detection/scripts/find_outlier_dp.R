# find outliers from linear regression between Mean and SD of region depth
# Input: pooled region depth
# Output: outliers individuals and regression graph

find.outlier.dp <- function(data, figure.file) {
    require(car, quietly=T)
    fit <- lm(SD ~ Mean, data=data)
    test.result <- outlierTest(fit, cutoff=0.01)
    outlier <- names(test.result$rstudent[test.result$rstudent > 0])
    png(figure.file)
    plot(data[,c("Mean","SD")])
    if (length(outlier) == 0) {
        message("\nNo samples have abnormal depth")
    } else {
        message("\nProblematic samples with abnormal mean and standard deviation of depth: ")
        data[outlier, "Bonferroni.p"] = test.result$bonf.p[outlier]
        print(data[outlier,c("Mean","SD","Bonferroni.p")])
        points(data[outlier,c("Mean","SD")],col="red")
    }
    abline(fit)
    dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
file.list <- args[1]
result.file <- args[2]
pngfile <- args[3]

require(data.table, quietly=T)
require(parallel, quietly=T)
files <- scan(file.list, character(), quiet=T)

# start parallel computing
cl <- makeCluster(detectCores())
region.mean <- rbindlist(parLapply(cl, files, fread))
region.mean <- region.mean[,2:(dim(region.mean)[2] - 1)]
result <- sapply(region.mean, function(x) {c(Mean = mean(x[x<100], na.rm=T), SD = sd(x[x<100], na.rm=T))})
result <- as.data.frame(t(result))
result <- result[order(rownames(result)), ]
#end parallel computing
stopCluster(cl)

find.outlier.dp(result, pngfile)
write.csv(result, result.file, quote=F)
