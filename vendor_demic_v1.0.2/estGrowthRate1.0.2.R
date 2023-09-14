library("demic")

args = commandArgs(trailingOnly=TRUE)

print(args[1])
X <- read.csv(args[1], header = FALSE, stringsAsFactors = TRUE)
colnames(X) <- c("logCov", "GC", "sample", "contig", "length")

output <- demic::estPTR(X, max_candidate_iter=args[3])

write.csv(output, file = paste(args[2], args[4], sep = "/"))