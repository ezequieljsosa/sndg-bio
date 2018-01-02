#!/bin/bash
#https://nelnet.org/docker/2017/03/23/Docker-Multiple-Commands-at-Run/
#do some stuff
#
env KEY_MODELLER=$1 dpkg -i /app/modeller_9.19-1_amd64.deb
bash -c "${@:2}"