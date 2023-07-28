#!/bin/bash

#Variables section
caseName=b10um

function makeCase(){
    mkdir -p ${caseName}

    cp -r src/0.orig    ${caseName}
    cp -r src/constant  ${caseName}
    cp -r src/system    ${caseName}
}

makeCase

echo -e "\nCreation done!\n"

exit 0
