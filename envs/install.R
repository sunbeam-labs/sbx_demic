if("lme4" %in% rownames(installed.packages()) == FALSE) {
    install.packages("lme4", repos="https://cloud.r-project.org")
}
if("factominer" %in% rownames(installed.packages()) == FALSE) {
    install.packages("factominer", repos="https://cloud.r-project.org")
}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {
    install.packages("reshape2", repos="https://cloud.r-project.org")
}