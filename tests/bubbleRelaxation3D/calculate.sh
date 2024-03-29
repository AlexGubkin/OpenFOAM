#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
caseName=case1

scaleSize=$(echo "0.000001" | bc -l)

cd ${caseName}

cp -r\
    0.orig\
    0

#Preparation mesh
foamDictionary  -entry decomposer -set scotch system/decomposeParDict
foamDictionary  -entry numberOfSubdomains -set 6 system/decomposeParDict
# foamDictionary  -entry simpleCoeffs/n -set "(1 2 3)" system/decomposeParDict
runApplication  blockMesh
runApplication  decomposePar
runParallel     renumberMesh -overwrite
runParallel     checkMesh -allGeometry -allTopology

#Calculation
runApplication  setFields
runParallel -a  setFields
runApplication  transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"
runParallel -a  transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"

foamDictionary  -entry endTime -set $(echo "0.005" | bc -l) system/controlDict
foamDictionary  -entry deltaT -set $(echo "0.00000001" | bc -l) system/controlDict
foamDictionary  -entry writeControl -set "adjustableRunTime" system/controlDict
# foamDictionary  -entry writeControl -set "timeStep" system/controlDict
foamDictionary  -entry writeInterval -set $(echo "0.000005" | bc -l) system/controlDict
# foamDictionary  -entry writeInterval -set 1 system/controlDict
foamDictionary  -entry adjustTimeStep -set "yes" system/controlDict
foamDictionary  -entry maxCo -set $(echo "0.2" | bc -l) system/controlDict
foamDictionary  -entry maxAlphaCo -set $(echo "0.2" | bc -l) system/controlDict

# runParallel     compressibleInterFoam
# runParallel     compressibleInterSSFFoam
runParallel     compressibleInterFSFFoam


exit 0
