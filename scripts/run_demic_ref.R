dir.create(snakemake@params[['output_dir']])

X <- read.csv(snakemake@input[['input']], stringsAsFactors = TRUE)
O <- demic::est_ptr(X)

write.table(O$all_ptr, snakemake@output[['all']], sep="\t", quote=FALSE, col.names=FALSE)
write.table(O$contigs_ptr, snakemake@output[['contig']], sep="\t", quote=FALSE, col.names=FALSE)
write.table(O$samples_ptr, snakemake@output[['sample']], sep="\t", quote=FALSE, col.names=FALSE)