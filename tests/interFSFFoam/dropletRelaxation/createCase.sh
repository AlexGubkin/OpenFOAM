#!/bin/bash

#Variables section
caseName=case1

function makeCase(){
    mkdir -p ${caseName}

    cp -r src/0.orig    ${caseName}
    cp -r src/constant  ${caseName}
    cp -r src/system    ${caseName}
    cp src/makeCloud.py ${caseName}
}

makeCase

echo -e "\nCreation done!\n"

exit 0
