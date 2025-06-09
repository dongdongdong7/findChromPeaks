chromData <- readRDS(file = "./inst/demo/chromData1.rds")
plot(x = chromData$rtime, y = chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime)
