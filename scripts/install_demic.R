#if ("demic" %in% rownames(installed.packages())) {
#    print("Already installed")
#} else {
    library(devtools)
    devtools::install_github("Ulthran/DEMIC")
#}

x <- data.frame()
write.table(x, file=snakemake@output[['out']], col.names=FALSE)