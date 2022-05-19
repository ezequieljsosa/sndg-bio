from ezequieljsosa/sndg-bio
RUN pip3 uninstall -y sndg-bio
RUN git clone https://github.com/ezequieljsosa/sndg-bio.git
RUN pip3 install mongoengine
ENV PYTHONPATH $PYTHONPATH:/sndg-bio
RUN apt update && apt install -y wget
RUN cd sndg-bio && git pull

# docker build --build-arg http_proxy=http://proxy.fcen.uba.ar:8080 --build-arg ftp_proxy=http://proxy.fcen.uba.ar:8080 --build-arg https_proxy=http://proxy.fcen.uba.ar:8080 -t ezequieljsosa/sndg-bio .

