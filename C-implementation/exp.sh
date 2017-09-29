#! /bin/bash
./bin/sim $1 $2 $3 $4 &&
FOLDER=$(ls ./results/ | tail -n 1) &&
# ls -hals ./results/$FOLDER/sim.bin &&
./bin/binToTxt ./results/$FOLDER/result.bin &&
./src/txtToPng/imGen.R $1 ./results/$FOLDER/result.txt &&
gnome-open ./results/$FOLDER/result.png
