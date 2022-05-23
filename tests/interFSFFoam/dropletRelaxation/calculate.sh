#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
caseName=case1

scaleSize=$(echo "0.001" | bc -l)

cd ${caseName}

cp -r\
    0.orig\
    0

cp\
    constant/src/fvModels\
    constant/fvModels

cp\
    constant/src/cloudProperties\
    constant/cloudProperties

python makeCloud.py

mv cloudPositions constant/

# cp\
#     constant/src/cloudPositions\
#     constant/cloudPositions

# cp\
#     constant/src/dynamicMeshDict\
#     constant/dynamicMeshDict

# cp\
#     constant/src/thermophysicalProperties.gas.0\
#     constant/thermophysicalProperties.gas
# 
# cp\
#     constant/src/thermophysicalProperties.liquid.0\
#     constant/thermophysicalProperties.liquid

# cp\
#     system/src/controlDict.0\
#     system/controlDict
# 
# cp\
#     system/src/fvSolution.0\
#     system/fvSolution

#Preparation mesh
foamDictionary  -entry numberOfSubdomains -set 6 system/decomposeParDict
foamDictionary  -entry method -set scotch system/decomposeParDict
runApplication  blockMesh
# runApplication  decomposePar -copyZero
runApplication  decomposePar

# runParallel     snappyHexMesh -overwrite
# runParallel     patchSummary
# runApplication  reconstructParMesh -constant
# runApplication  reconstructPar -withZero
# rm log.decomposePar
# runApplication  decomposePar -force
runParallel     renumberMesh -overwrite

# for d in processor* ; do
# #     echo "${d}"
#     rm ${d}/constant/polyMesh/sets/protectedCells*
# done

# runParallel     createPatch -overwrite

# runApplication  extrudeMesh
# # rm log.decomposePar
# # runApplication  decomposePar -force
# # runParallel     redistributePar
# # runParallel     renumberMesh -overwrite -noFields
# runParallel     patchSummary
# runParallel     checkMesh -allGeometry -allTopology

#Calculation
runParallel     setFields
runApplication  transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"
rm log.transformPoints
runParallel     transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"

foamDictionary  -entry endTime -set $(echo "10.0" | bc -l) system/controlDict
foamDictionary  -entry writeInterval -set $(echo "0.001" | bc -l) system/controlDict

runParallel     compressibleInterFSFFoam

# cp\
#     constant/src/dynamicMeshDict\
#     constant/dynamicMeshDict
# 
# rm log.thermoCapillarityInterFoam
# 
# for d in processor* ; do
#     foamDictionary  -entry deltaT -set $(echo "0.0000000001" | bc -l) ${d}/$(foamListTimes -processor -latestTime)/uniform/time
#     foamDictionary  -entry deltaT0 -set $(echo "0.0000000001" | bc -l) ${d}/$(foamListTimes -processor -latestTime)/uniform/time
# done
# 
# foamDictionary  -entry endTime -set $(echo "10.0" | bc -l) system/controlDict
# foamDictionary  -entry writeInterval -set $(echo "0.001" | bc -l) system/controlDict
# 
# runParallel     thermoCapillarityInterFoam

exit 0
