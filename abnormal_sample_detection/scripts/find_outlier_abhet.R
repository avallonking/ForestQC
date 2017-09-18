find.outlier.abhet <- function(data, upperbound=0.6, lowerbound=0.4) {
    outlier.idx <- which(data$ABHet >= upperbound | data$ABHet <= lowerbound)
    if (length(outlier.idx) == 0) {
        message("No samples with abnormal ABHet")
    } else {
        message("\nProblematic samples with abnormal ABHet: ")
        print(data[outlier.idx, "ABHet", drop=F])
    }
}

args <- commandArgs(trailingOnly = TRUE)
file.list <- args[1]
outfile <- args[2]
files <- scan(file.list, character(), quiet=T)

require(data.table, quietly=T)
#require(parallel)
#cl <- makeCluster(detectCores())
#abhet <- do.call('+', lapply(files, function(x) data.frame(fread(x), row.names=1, check.names=F)))
#stopCluster()
for (file in files) {
    temp <- data.frame(fread(file), row.names=1, check.names=F)
    if (!exists("abhet")) {
        abhet <- temp
    } else {
        abhet <- abhet + temp
    }
}

abhet <- as.data.frame(t(abhet))
abhet <- abhet[order(rownames(abhet)), ]
abhet$ABHet <- abhet$REF / (abhet$REF + abhet$ALT)
find.outlier.abhet(abhet, upperbound = 0.6, lowerbound=0.4)
write.csv(abhet[, "ABHet", drop=F], outfile, quote=F, row.names=T)
