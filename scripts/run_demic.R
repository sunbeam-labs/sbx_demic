files <- list.files(path=snakemake@input[['input']], pattern="*.cov3", full.names=TRUE, recursive=FALSE)
dir.create(snakemake@output[['out']])
lapply(files, function(x) {
    print(paste("Working on: ", x))
    tryCatch(
        {
            X <- read.csv(x, stringsAsFactors = TRUE)
            O <- demic::est_ptr(X)
            
            write.table(O$all_ptr, paste(snakemake@output[['out']], "/", tools::file_path_sans_ext(basename(x)), ".all.ptr", sep=""), sep="\t", quote=FALSE, col.names=FALSE)
            write.table(O$contigs_ptr, paste(snakemake@output[['out']], "/", tools::file_path_sans_ext(basename(x)), ".contig.ptr", sep=""), sep="\t", quote=FALSE, col.names=FALSE)
            write.table(O$samples_ptr, paste(snakemake@output[['out']], "/", tools::file_path_sans_ext(basename(x)), ".sample.ptr", sep=""), sep="\t", quote=FALSE, col.names=FALSE)
        },
        error=function(e) {
            print(e)
        },
        warning=function(w) {
            print(w)
            return(NA)
        }
    )
})
