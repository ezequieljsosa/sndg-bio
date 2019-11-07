FROM pescobar/openstructure
#https://hub.docker.com/r/pescobar/openstructure/dockerfile
#https://git.scicore.unibas.ch/schwede/QMEAN/
WORKDIR /opt
RUN wget https://git.scicore.unibas.ch/schwede/QMEAN/-/archive/master/QMEAN-master.tar.gz 
RUN tar xfv  QMEAN-master.tar.gz && mv QMEAN-master qmean && rm /opt/QMEAN-master.tar.gz
RUN cd qmean && mkdir build &&  cd build && cmake ../ -DOST_ROOT=/usr/local/lib64/ -DOPTIMIZE=1 -DOST_INCLUDE_DIR=/usr/local/include/
RUN cd qmean/build && make 
RUN yum install -y python-matplotlib
PYTHONPATH=/opt/qmean/build/stage/lib64/python2.7/site-packages/
ENTRYPOINT bash
