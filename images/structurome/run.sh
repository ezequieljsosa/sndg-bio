#!/bin/bash
#https://nelnet.org/docker/2017/03/23/Docker-Multiple-Commands-at-Run/
#https://stackoverflow.com/questions/16988427/calling-one-bash-script-from-another-script-passing-it-arguments-with-quotes-and

if [ "$1" == "xxx" ]
then
    echo "A valid Modeller key must be entered as first parameter"
    exit 1
fi

if [ "$1" == "createdb.py" ]
then
    /app/createdb.py
    exit 0
fi

if [ "$1" == "x" ]
then
    theargs="'$2'"
    for i in "${@:3}" ; do
       theargs="${theargs} '$i'"
    done

    bash -c "$theargs"
    exit 0
fi


env KEY_MODELLER=$1 dpkg -i /app/modeller_9.19-1_amd64.deb
theargs="'$2'"
for i in "${@:3}" ; do
   theargs="${theargs} '$i'"
done

bash -c "$theargs"