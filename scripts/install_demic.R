#if ("demic" %in% rownames(installed.packages())) {
#    writeLines(c("Already installed"), snakemake@log[['log']])
#} else {
    writeLines(c("Installing demic"), snakemake@log[['log']])
    #install.packages("rlang", repos = "http://cran.us.r-project.org")
    #install.packages("usethis", repos = "http://cran.us.r-project.org", dependencies=TRUE)
    #install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies=TRUE)
    library(remotes)
    remotes::install_github("Ulthran/DEMIC@15-release-demic-100")
#}

library(demic)
x <- data.frame()
write.table(x, file=snakemake@output[['out']], col.names=FALSE)