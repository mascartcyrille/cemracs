#!/usr/bin/Rscript
rm(list=ls())
args=commandArgs()
if(interactive()) {
	args<-c("1","2","3","4","5","result.txt")
}
for( i in 6:length(args) ) {
	timesFile=args[i]
	# pictureFile=stringr::str_replace(stringr::str_replace(stringr::str_replace(timesFile, ".csv", ".png"), "spikeTimes", "spikes"), "/Times/", "/Pictures/")
	pictureFile = "result.png"
	Times=read.csv(args[i], header=F)[[1]]
	S=1:length(Times)
	png(pictureFile)
	plot(Times,S)
	dev.off()
}