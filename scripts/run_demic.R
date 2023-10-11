files <- list.files(path=snakemake@input[['input']], pattern="*.cov3", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
    X <- read.csv(x, header = FALSE, stringsAsFactors = TRUE)
    colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
    O <- demic::est_ptr(X)
    
    dir.create(snakemake@output[['out']])
    write.table(O, paste(snakemake@output[['out']], "/", tools::file_path_sans_ext(basename(x)), ".ptr", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
})