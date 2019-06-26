FROM ubuntu:16.04

WORKDIR /app
RUN  apt-get update && apt-get install -y  build-essential curl file git  software-properties-common
#RUN add-apt-repository ppa:jonathonf/python-3.6
RUN apt-get update && apt-get install -y python3-dev python3-pip python-dev python-pip
RUN useradd --system -s /sbin/nologin someuser && mkdir /home/someuser && chown someuser:someuser  /home/someuser
USER someuser
ENV PATH=/root/.linuxbrew/Homebrew/Library/Homebrew/vendor/portable-ruby/current/bin:$PATH
RUN USER=someuser bash -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
USER root
RUN pip3 install numpy
RUN echo 'eval $(/home/someuser/.linuxbrew/bin/brew shellenv)' >>~/.profile
RUN eval $(/home/someuser/.linuxbrew/bin/brew shellenv) && brew install bwa samtools bcftools blast
#RUN pip3 install --upgrade pip
RUN apt install  -y muscle
#ENV C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/python3.6m/
RUN apt install -y libz-dev  zlib1g-dev libssl-dev libcurl4-openssl-dev bzip2 liblzma-dev libbz2-dev


RUN apt-get install -y python3-lxml python3-six




RUN apt install -y python3-pyqt5
RUN apt install -y python3-pyqt5.qtopengl
RUN apt install -y python3-pyqt5.qtsvg

RUN pip3 install ete3

ENV DISPLAY :0
RUN apt install -y xvfb
#RUN apt-get install -qqy x11-apps nano
RUN mkdir /usr/lib/x86_64-linux-gnu/fonts/
RUN apt -y install ttf-liberation
RUN ln  /usr/share/fonts/truetype/liberation/Liberation*  /usr/lib/x86_64-linux-gnu/fonts/

RUN pip2 install numpy
RUN pip2 install matplotlib==2 scipy==1.2
RUN pip2 install pyvcf
RUN pip3 install pyvcf
RUN pip2 install VCF-kit
RUN pip3 install goatools
RUN pip3 install sndg-bio==0.1.26
RUN sed -i 's/np.set_printoptions(threshold=np.nan)/#np.set_printoptions(threshold=np.nan)/' /usr/local/lib/python2.7/dist-packages/vcfkit/utils/vcf.py

#ENTRYPOINT ["/app/run.sh"]
WORKDIR /out
CMD ["bash"]
