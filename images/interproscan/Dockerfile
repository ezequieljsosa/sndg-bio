FROM ubuntu:18.04

RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y libdw-dev  build-essential python python3 openjdk-11-jdk && \
        apt-get clean && apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
