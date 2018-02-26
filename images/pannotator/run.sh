#!/bin/bash
#https://nelnet.org/docker/2017/03/23/Docker-Multiple-Commands-at-Run/
#https://stackoverflow.com/questions/16988427/calling-one-bash-script-from-another-script-passing-it-arguments-with-quotes-and

if [ "$1" == "createdb.py" ]
then
    /app/p_procariota/createdb.py
    exit 0
fi

theargs="'$1'"
for i in "${@:2}" ; do
   theargs="${theargs} '$i'"
done

bash -c "$theargs"