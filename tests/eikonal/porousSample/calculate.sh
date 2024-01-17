#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
scaleForFolderNames=4
scaleForCalculations=12

scaleSize=$(echo "scale=${scaleForCalculations}; 0.001/1.0;" | bc -l)

WD=$(pwd)
srcDir=${WD}/src
solutionDir=${WD}
backgroundMeshDir=${WD}/backgroundMesh

# export TMP=${WD}
# export TEMP=${WD}
# export TMPDIR=${WD}
#
# mktemp -d

caseDir=${solutionDir}/case1

mkdir -p ${caseDir}

cd ${caseDir}

echo -e "\n######################################################################"
echo -e "Case:" ${caseDir}
echo -e "######################################################################\n"

#Copy sources to case directory
cp -ru ${srcDir}/0/ ${caseDir}
cp -ru ${srcDir}/system/ ${caseDir}

#Copy background mesh to case directory
cp -ru ${backgroundMeshDir}/constant/polyMesh ${caseDir}/constant
cp -ru ${backgroundMeshDir}/processor* ${caseDir}
for d in processor* ; do
    rm ${d}/constant/polyMesh/*Level*
done
runApplication -o\
    foamDictionary\
        -entry entry0/zPlane/type\
        -set empty\
        constant/polyMesh/boundary
runApplication -a\
    foamDictionary\
        -entry entry0/dzPlane/type\
        -set empty\
        constant/polyMesh/boundary

#Preparation STL-files for snappyHexMesh
mkdir -p ${caseDir}/constant/triSurface/
tar -xzf ${srcDir}/skeletonSample3D.tar.gz -C ${caseDir}
ln -sf ${caseDir}/openPoresSkeleton3D.stl ${caseDir}/constant/triSurface/rockSkeletonWalls.stl

#             runApplication -o\
#                 surfaceFeatures\
#                     -case ${caseDir}\
#                     -dict ${caseDir}/system/snappyHexMeshSrc/surfaceFeaturesDict
#             runApplication -o\
#                 surfaceFeatureConvert\
#                     -case ${caseDir}\
#                     ${caseDir}/constant/triSurface/rockSkeletonWalls.eMesh\
#                     ${caseDir}/constant/triSurface/rockSkeletonWalls.vtk

# Run snappyHexMesh
runApplication -a\
    foamDictionary\
        -entry numberOfSubdomains\
        -set 60\
        system/decomposeParDict
runApplication -a\
    foamDictionary\
        -entry method\
        -set scotch\
        system/decomposeParDict
runParallel -o\
    snappyHexMesh -overwrite\
        -case ${caseDir}\
        -dict ${caseDir}/system/snappyHexMeshSrc/snappyHexMeshDict.porousSample
runParallel -o\
    transformPoints\
        -case ${caseDir}\
        "scale=(${scaleSize} ${scaleSize} ${scaleSize})"
runParallel -o\
    checkMesh -allTopology -allGeometry\
        -case ${caseDir}
runParallel -a\
    foamDictionary\
        -entry entry0/zPlane/type\
        -set empty\
        constant/polyMesh/boundary
runParallel -a\
    foamDictionary\
        -entry entry0/dzPlane/type\
        -set empty\
        constant/polyMesh/boundary
for d in processor* ; do
    cp -ru 0 ${d}
done
runParallel -o\
    patchSummary\
        -case ${caseDir}

#             runApplication -o\
#                 reconstructParMesh\
#                     -constant\
#                     -case ${caseDir}

runParallel -o\
    HeleShawSimpleFoam\
        -case ${caseDir}

echo -e "\nCalculation done!\n"

exit 0
