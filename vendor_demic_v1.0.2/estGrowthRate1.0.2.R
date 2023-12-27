library("demic")

args = commandArgs(trailingOnly=TRUE)

print(args[1])
X <- read.csv(args[1], stringsAsFactors = TRUE)
colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")

output <- demic::est_ptr(X, max_candidate_iter=args[3])

write.csv(output, file = paste(args[2], args[4], sep = "/"))