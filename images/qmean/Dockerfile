#FROM registry.scicore.unibas.ch/schwede/openstructure:2.0.0
FROM registry.scicore.unibas.ch/schwede/openstructure
#https://git.scicore.unibas.ch/schwede/openstructure/tree/2.0.0/docker


#https://git.scicore.unibas.ch/schwede/QMEAN/
WORKDIR /opt
RUN wget https://git.scicore.unibas.ch/schwede/QMEAN/-/archive/master/QMEAN-master.tar.gz 
RUN tar xfv  QMEAN-master.tar.gz && mv QMEAN-master qmean && rm /opt/QMEAN-master.tar.gz
#RUN yum install -y cmake3
#ENV OST_ROOT=/usr/local/src/openstructure-1.8.0-source
#WORKDIR /usr/local/src
#RUN wget -O openstructure-1.8.0.tar.gz "https://git.scicore.unibas.ch/schwede/openstructure/-/archive/1.8.0/openstructure-1.8.0.tar.gz"
#RUN tar -xf openstructure-1.8.0.tar.gz --strip-components=1 -C /usr/local/src/openstructure-1.8.0-source
#WORKDIR /usr/local/src/openstructure-${OPENSTRUCTURE_VERSION}-source/build
#RUN cmake3 .. \
#   -DENABLE_MM=1 \
#   -DOPEN_MM_LIBRARY=/usr/local/openmm/lib/libOpenMM.so \
#   -DOPEN_MM_PLUGIN_DIR=/usr/local/openmm/lib/plugins \
#   -DOPEN_MM_INCLUDE_DIR=/usr/local/openmm/include \
#   -DCOMPILE_TMTOOLS=1 \
#   -DENABLE_GFX=ON \
#   -DENABLE_GUI=ON \
#   -DUSE_NUMPY=1 \
#   -DUSE_RPATH=1 \
#   -DEIGEN3_INCLUDE_DIR=/usr/local/include/Eigen/ \
#   -DFFTW_LIBRARY=/usr/lib64/libfftw3f.a \
#   -DQT_QMAKE_EXECUTABLE=/usr/lib64/qt4/bin/qmake \
#   -DOPTIMIZE=1 && \
#   make -j ${CPUS_FOR_MAKE} && \
#   make check && \
#   make install

WORKDIR /opt

#ENV PYTHONPATH=/opt/qmean/build/stage/lib64/python3.6/site-packages/

#ENV OST_ROOT=/usr/local/src/openstructure-1.8.0-source/build/stage/

#RUN yum install -y python-devel
#RUN yum install -y python3-devel

RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.1/cmake-3.18.1.tar.gz && tar xfv cmake-3.18.1.tar.gz
RUN apt update && apt install -y libssl-dev
RUN cd cmake-3.18.1 && ./configure  && make
RUN cd cmake-3.18.1 && make install


ENV PYTHONPATH=$PYTHONPATH:/usr/include/x86_64-linux-gnu/python3.6m/
RUN cd qmean && mkdir build &&  cd build && cmake ../ -DOST_ROOT=/usr/local/ \
    -DOPTIMIZE=1 -DOST_INCLUDE_DIR=/usr/local/src/openstructure-1.8.0-source/build/stage/include \
    -DPYTHON_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython3.6m.so \
    -DPYTHON_INCLUDE_PATH=/usr/include/python3.6m/



RUN cd qmean/build && make
RUN apt install -y python3-matplotlib
ENV PYTHONPATH=$PYTHONPATH:/opt/qmean/build/stage/lib64/python3.6/site-packages/
#ENV PYTHONPATH=$PYTHONPATH:/usr/local/src/openstructure-1.8.0-source/build/stage/lib64/python2.7/site-packages/:/usr/local/src/openstructure-1.8.0-source/build/stage/
#ENV PYTHONPATH=$PYTHONPATH:/usr/local/src/openstructure-1.8.0-source/build/stage/lib64/python2.7/site-packages/:/usr/local/src/openstructure-1.8.0-source/build/stage/
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/qmean/build/stage/lib64/
COPY qmean /usr/local/bin
RUN chmod +x /usr/local/bin/qmean
ENTRYPOINT []
CMD qmean
