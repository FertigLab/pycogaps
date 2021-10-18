# library(CoGAPS)

R.utils::sourceDirectory("/Users/ashleytsang/pycogaps/src/CoGAPS/src")
R.utils::sourceDirectory("/Users/ashleytsang/pycogaps/src/CoGAPS/R")


path <- "/Users/ashleytsang/pycogaps/data/GSE98638_HCC.TCell.S5063.count.txt"
counts <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
counts <- counts[-c(1, 2)]


params <- new("CogapsParams")
params <- setParam(params, "nPatterns", 3)
params <- setParam(params, "sparseOptimization", TRUE)
params <- setParam(params, "nIterations", 10000)
params <- setParam(params, "seed", 0)



start_time <- Sys.time()
CoGAPS(counts, params)
end_time <- Sys.time()

diff <- end_time - start_time
diff