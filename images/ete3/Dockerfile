FROM ubuntu:16.04
RUN apt update && apt install -y python3-pip wget curl unzip htop
RUN apt install -y python3-pip wget curl unzip htop
RUN apt install -y libxml2-dev libxslt-dev
RUN apt install -y python3-pyqt5 python3-pyqt5.qtopengl python3-pyqt5.qtsvg
RUN apt install -y python3-six python3-numpy
RUN ln -snf /usr/share/fonts/truetype/dejavu/ /usr/lib/x86_64-linux-gnu/fonts
RUN pip3 install ete3
RUN apt install -y python3-pandas

#COPY setup.py /sndg-bio/setup.py
#WORKDIR /sndg-bio
#RUN pip3 install . && pip uninstall -y `python3 setup.py --name`
#COPY . /sndg-bio/

COPY render_tree.py /usr/local/bin/render_tree.py

#ENV PYTHONPATH /sndg-bio/
ENV QT_QPA_PLATFORM offscreen

