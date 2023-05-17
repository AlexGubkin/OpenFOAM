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
runApplication -a\
    foamDictionary  -entry decomposer -set scotch system/decomposeParDict
runApplication -a\
    foamDictionary  -entry numberOfSubdomains -set 40 system/decomposeParDict
# runApplication -a\
#     foamDictionary  -entry simpleCoeffs/n -set "(7 7 1)" system/decomposeParDict

runApplication  blockMesh
runApplication  decomposePar -copyZero
runParallel     snappyHexMesh -overwrite
# runParallel     renumberMesh -overwrite
runParallel     extrudeMesh
runParallel     createPatch -overwrite
runParallel     checkMesh -allGeometry -allTopology

runParallel -a\
    foamDictionary  constant/polyMesh/boundary -entry entry0/back/type -set empty
runParallel -a\
    foamDictionary  constant/polyMesh/boundary -entry entry0/front/type -set empty
runApplication -a\
    foamDictionary  constant/polyMesh/boundary -entry entry0/back/type -set empty
runApplication -a\
    foamDictionary  constant/polyMesh/boundary -entry entry0/front/type -set empty

runApplication reconstructParMesh -constant -withZero

#Calculation
runApplication  setFields
runParallel -a  setFields
# runApplication  transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"
# runParallel -a  transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"

foamDictionary  -entry endTime -set $(echo "20.0" | bc -l) system/controlDict
foamDictionary  -entry deltaT -set $(echo "0.00000001" | bc -l) system/controlDict
foamDictionary  -entry writeControl -set "adjustableRunTime" system/controlDict
# foamDictionary  -entry writeControl -set "timeStep" system/controlDict
foamDictionary  -entry writeInterval -set $(echo "0.01" | bc -l) system/controlDict
# foamDictionary  -entry writeInterval -set 1 system/controlDict
foamDictionary  -entry adjustTimeStep -set "yes" system/controlDict
foamDictionary  -entry maxCo -set $(echo "0.2" | bc -l) system/controlDict
foamDictionary  -entry maxAlphaCo -set $(echo "0.2" | bc -l) system/controlDict

runParallel     interHeleShawFoam


exit 0
