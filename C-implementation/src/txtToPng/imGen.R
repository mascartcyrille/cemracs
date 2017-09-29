#!/usr/bin/Rscript
rm( list = ls() )
args = commandArgs()
if( interactive() ) {
	args<-c( "1", "2", "3", "4", "5", "10000", "result.txt" )
}
for( i in 7:length( args ) ) {
	timesFile	=	args[i]
	pictureFile	=	stringr::str_replace(timesFile, ".txt", ".png")
	Times		=	read.csv(args[i], header=F)[[1]]
	#S			=	(1:length( Times )) / length( Times )
	S			=	1:length( Times )
	
	png( pictureFile )
	plot( Times, S )
	dev.off()
}
