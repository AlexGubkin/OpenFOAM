#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
caseName=case1

scaleSize=$(echo "0.001" | bc -l)

cd ${caseName}

cp -r\
    0.orig\
    0

#Preparation mesh
foamDictionary  -entry numberOfSubdomains -set 8 system/decomposeParDict
foamDictionary  -entry method -set scotch system/decomposeParDict
runApplication  blockMesh
runApplication  decomposePar
runParallel     renumberMesh -overwrite
runParallel     checkMesh -allGeometry -allTopology

#Calculation
runParallel     setFields
runApplication  transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"
rm log.transformPoints
runParallel     transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"

foamDictionary  -entry endTime -set $(echo "0.1" | bc -l) system/controlDict
foamDictionary  -entry writeInterval -set $(echo "0.001" | bc -l) system/controlDict

runParallel     interFSFFoam

exit 0
