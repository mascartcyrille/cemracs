#!/bin/bash
while test $# -gt 0; do
	case "$1" in
		-c|--compile)
			mkdir -p build/
			javac -d ./build/ -cp ./src/:./jars/ssj-2.5.jar src/TestGUI/*.java src/TestSingleModel/*.java src/Util/*.java
			cd build/
			jar -cfe ../jars/cl.jar TestSingleModel.Simulator TestSingleModel/*.class Util/*.class ../jars/ssj-2.5.jar
			jar -cfe ../jars/g.jar TestGUI.execGUI TestGUI/*.class TestSingleModel/*.class Util/*.class ../jars/ssj-2.5.jar
			cd ..
			shift
			;;
		-h|--help)
			echo "sim -h|--help"
			echo "sim -c|--compile"
			echo "sim [-c|--compile] -e|--exec [-g] neuronNumber connectionProbability iterationNumber"
			exit 0
			;;
		-e|--exec)
			shift
			exec="cl.jar TestSingleModel.Simulator"
			if test $# -gt 2; then
				if [ $1 = "-g" ]
				then
					shift
					exec="g.jar TestGUI.execGUI"
				fi
				java -classpath jars/ssj-2.5.jar:jars/$exec $1 $2 $3
				exit 0
			else
				echo "Bad arguments"
				exit 1
			fi
			;;
		*)
			echo "Bad arguments"
			exit 1
			;;
	esac
done
