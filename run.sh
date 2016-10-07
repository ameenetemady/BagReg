#!/bin/bash
bClassPath=./lib/weka.jar
classPath=./bin:$bClassPath

#compile
#javac -sourcepath src/proteinidentification -classpath $bClassPath  -d bin src/proteinidentification/*.java

pep_identification=importdata/Sigma_49.txt
prot_res=importdata/Sigma_49_result.csv
prot_refListFile=ref/Sigma_49_reference.csv

java -cp $classPath  proteinidentification.RAAlgorithm $pep_identification $prot_res
python3 getAUC.py $prot_refListFile $prot_res 
