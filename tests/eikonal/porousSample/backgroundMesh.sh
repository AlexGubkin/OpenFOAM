#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

WD=$(pwd)
srcDir=${WD}/src
backgroundMeshDir=${WD}/backgroundMesh

mkdir -p ${backgroundMeshDir}/system

cd ${backgroundMeshDir}

cp -u\
    ${srcDir}/system/blockMeshDict\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/snappyHexMeshSrc/snappyHexMeshDict.backgroundMesh\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/snappyHexMeshSrc/locationInMesh\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/controlDict\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/decomposeParDict\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/extrudeMeshDict\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/fvSchemes\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/fvSolution\
    ${backgroundMeshDir}/system
cp -u\
    ${srcDir}/system/ROI\
    ${backgroundMeshDir}/system

#Preparation background mesh for snappyHexMesh
runApplication -o\
    blockMesh\
        -case ${backgroundMeshDir}\

#Run snappyHexMesh
runApplication -o\
    foamDictionary\
        -entry numberOfSubdomains\
        -set 60 system/decomposeParDict
runApplication -a\
    foamDictionary\
        -entry method\
        -set scotch system/decomposeParDict
runApplication -o\
    decomposePar\
        -force\
        -case ${backgroundMeshDir}
runParallel -o\
    snappyHexMesh\
        -overwrite\
        -case ${backgroundMeshDir}\
        -dict ${backgroundMeshDir}/system/snappyHexMeshDict.backgroundMesh
runParallel -o\
    extrudeMesh\
        -case ${backgroundMeshDir}
runParallel -o\
    checkMesh\
        -allTopology\
        -allGeometry\
        -case ${backgroundMeshDir}
runParallel -a\
    foamDictionary\
        -entry entry0/zPlane/type\
        -set patch constant/polyMesh/boundary
runParallel -a\
    foamDictionary\
        -entry entry0/dzPlane/type\
        -set patch constant/polyMesh/boundary
runParallel -o\
    patchSummary\
        -case ${backgroundMeshDir}
