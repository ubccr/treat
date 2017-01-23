#!/bin/bash

TREAT_DIR='./.treat-release'
VERSION=`git describe --long --tags --dirty --always | sed -e 's/^v//'`

for os in linux darwin windows
do
    for arch in amd64 386
    do
        NAME=treat-${VERSION}-${os}-${arch}
        REL_DIR=${TREAT_DIR}/${NAME}
        cd ./cmd/treat && GOOS=$os GOARCH=$arch go build -ldflags "-X main.TreatVersion=$VERSION" .
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
done
