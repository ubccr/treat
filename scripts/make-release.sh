#!/bin/bash

source $GOROOT/src/golang-crosscompile/crosscompile.bash
TREAT_DIR='./.treat-release'
VERSION=`grep Version cmd/treat/main.go | egrep -o '[0-9]\.[0-9]\.[0-9]'`

for arch in linux-amd64 darwin-amd64 windows-amd64 windows-386
do
    NAME=treat-${VERSION}-${arch}
    REL_DIR=${TREAT_DIR}/${NAME}
    GO_BIN="go-${arch}"
    cd ./cmd/treat && ${GO_BIN} build . 
    cd ../../
    rm -Rf ${TREAT_DIR}
    mkdir -p ${REL_DIR}
    cp ./cmd/treat/treat* ${REL_DIR}/ 
    cp ./README.rst ${REL_DIR}/ 
    cp ./AUTHORS.rst ${REL_DIR}/ 
    cp ./ChangeLog.rst ${REL_DIR}/ 
    cp ./LICENSE ${REL_DIR}/ 
    cp -R ./docs ${REL_DIR}/ 
    cp -R ./examples ${REL_DIR}/ 
    cp -R ./cmd/treat/templates ${REL_DIR}/ 

    cd ${TREAT_DIR} && zip -r ${NAME}.zip ${NAME}
    mv  ${NAME}.zip ../
    cd ../
    rm -Rf ${TREAT_DIR}
    rm ./cmd/treat/treat*
done
