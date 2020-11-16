FROM ubuntu:16.04

WORKDIR /app


# docker build -t structurome --build-arg=http_proxy=http://proxy.fcen.uba.ar:8080 --build-arg=https_proxy=http://proxy.fcen.uba.ar:8080  --build-arg=ftp_proxy=http://proxy.fcen.uba.ar:8080 .

RUN apt-get update && apt-get install -y libnetcdf-dev \
    libfftw3-dev libtiff5-dev libeigen3-dev wget

RUN wget https://netix.dl.sourceforge.net/project/fpocket/fpocket2.tar.gz -O fpocket2.tar.gz && tar xfv fpocket2.tar.gz \
    && rm fpocket2.tar.gz &&  cd /app/fpocket2 && sed -i 's|$(LINKER) $(LFLAGS) $^ -o $@|$(LINKER) $^ -o $@ $(LFLAGS)|' makefile && make && make install



WORKDIR /out

CMD ["bash"]
