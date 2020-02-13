args = commandArgs(trailingOnly=TRUE) 
IP <- installed.packages()
write(c(R.version$major, R.version$minor, IP["lme4","Version"], IP["FactoMineR","Version"]), args[1])