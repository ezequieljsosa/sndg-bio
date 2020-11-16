
docker build -t "ezequieljsosa/reads_processing"  --build-arg http_proxy="http://proxy.fcen.uba.ar:8080" \
  --build-arg https_proxy="http://proxy.fcen.uba.ar:8080" \
  --build-arg ftp_proxy="http://proxy.fcen.uba.ar:8080" .