#install.packages("strategusMultiverse_0.1.2.tar.gz", repos = NULL, type = "source")
library(strategusMultiverse)

resultsPath <- "../resultsFolder/results_folder"
mv <- readStrategusMultiverse(resultsPath)
mv


multiverse(resultsPath, type = "specification_curve")
multiverse(resultsPath, type = "multiverse")
multiverse(resultsPath, type = "density")
multiverse(resultsPath, type = "volcano")
multiverse(resultsPath, type = "influence")
multiverse(resultsPath, type = "gate")
multiverse(resultsPath, type = "shiny")

# specification_curve, multiverse, density, volcano, influence, gate, shiny