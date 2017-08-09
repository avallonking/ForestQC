find.outlier.abhet <- function(data, upperbound=0.6, lowerbound=0.4) {
    outlier.idx <- which(data$ABHet >= upperbound | data$ABHet <= lowerbound)
    message("Problematic samples with abnormal ABHet: ")
    print(data[outlier.idx, "ABHet", drop=F])
}

args <- commandArgs(trailingOnly = TRUE)
file.list <- args[1]
outfile <- args[2]
files <- scan(file.list, character())
for (file in files) {
    temp <- read.csv(file, row.names=1, check.names=F)
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
abhet <- abhet["ABHet"]
write.table(abhet, outfile, quote=F, row.names=F)
