#if ("demic" %in% rownames(installed.packages())) {
#    writeLines(c("Already installed"), snakemake@log[['log']])
#} else {
    writeLines(c("Installing demic"), snakemake@log[['log']])
    library(remotes)
    remotes::install_github("Ulthran/DEMIC@master")
#}

library(demic)
x <- data.frame()
write.table(x, file=snakemake@output[['out']], col.names=FALSE)