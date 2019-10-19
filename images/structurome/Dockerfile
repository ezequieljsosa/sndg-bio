FROM ubuntu:16.04

WORKDIR /app


# docker build -t structurome --build-arg=http_proxy=http://proxy.fcen.uba.ar:8080 --build-arg=https_proxy=http://proxy.fcen.uba.ar:8080  --build-arg=ftp_proxy=http://proxy.fcen.uba.ar:8080 .

RUN apt-get update && apt-get install -y python-tk python-pip libnetcdf-dev \
   python-biopython cmake libgtk2.0-dev ncbi-blast+ libboost-all-dev libfftw3-dev libtiff5-dev libeigen3-dev wget


RUN pip install numpy scipy sklearn matplotlib tqdm


ADD /data/ost /app/ost
ADD /data/qmean /app/qmean
ADD /data/dssp /app/
ADD /data/modeller_9.19-1_amd64.deb /app/




RUN cd /app/ost && cmake . -DOPTIMIZE=ON -DENABLE_INFO=OFF -DPYTHON_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython2.7.so \
    && make && mv stage /opt/ost
RUN cd /app/qmean && cmake . -DOST_ROOT=/opt/ost  -DOPTIMIZE=1 -DPYTHON_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython2.7.so \
    && make && mv stage /opt/qmean

ADD /data/fpocket-3.1 /app/fpocket-3.1
RUN rm /app/fpocket-3.1/src/mdpocket.c
RUN cd /app/fpocket-3.1 && make && make install

ADD run.sh /app
ADD createdb.py /app
RUN chmod +x /app/run.sh && ln -s /app/dssp /usr/local/bin/dssp && rm -r /app/ost /app/qmean

RUN pip install sndg-bio==0.1.16


ENV PYTHONPATH=/app/structurome/:/opt/qmean/lib64/qmean:/opt/ost/lib64/python2.7/site-packages/:/opt/ost/
ENV LD_LIBRARY_PATH=/opt/ost/lib64/:/opt/qmean/lib64/

WORKDIR /out

ENTRYPOINT ["/app/run.sh"]
CMD ["xxx"]
