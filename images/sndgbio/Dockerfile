from ezequieljsosa/sndg-bio
RUN pip3 uninstall -y sndg-bio
RUN git clone https://github.com/ezequieljsosa/sndg-bio.git
RUN pip3 install mongoengine
ENV PYTHONPATH $PYTHONPATH:/sndg-bio
RUN apt update && apt install -y wget
