#!/bin/bash
while test $# -gt 0; do
	case "$1" in
		-c|--compile)
			mkdir -p build/
			javac -d ./build/ -cp ./src/:./jars/ssj-2.5.jar src/TestSingleModel/*.java src/Util/*.java
			cd build/
			jar -cfe ../jars/cl.jar TestSingleModel.Simulator TestSingleModel/*.class Util/*.class ../jars/ssj-2.5.jar
			cd ..
			shift
			;;
		-h|--help)
			echo "sim -h|--help"
			echo "sim -c|--compile"
			echo "sim [-c|--compile] -e|--exec neuronNumber connectionProbability iterationNumber"
			exit 0
			;;
		-e|--exec)
			shift
			exec="cl.jar TestSingleModel.Simulator"
			if test $# -gt 2; then
				java -classpath jars/ssj-2.5.jar:jars/$exec $1 $2 $3
				exit 0
			else
				echo "Bad arguments"
				exit 1
			fi
			;;
		*)
			echo "Bad arguments"
			echo "Try ./sim.sh [-h|--help]"
			exit 1
			;;
	esac
done
