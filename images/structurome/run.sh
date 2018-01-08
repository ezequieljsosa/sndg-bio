#!/bin/bash
#https://nelnet.org/docker/2017/03/23/Docker-Multiple-Commands-at-Run/

if [ "$1" == "xxx" ]
then
    echo "A valid Modeller key must be entered as first parameter"
    return 1
fi

if [ "$1" == "createdb.py" ]
then
    /app/createdb.py
    return 0
fi


env KEY_MODELLER=$1 dpkg -i /app/modeller_9.19-1_amd64.deb
bash -c "${@:2}"