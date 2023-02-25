#FROM bzhanglab/jvarkit
FROM jfx319/jvarkit:bam2raster

WORKDIR /git/jvarkit

RUN make fixvcfmissinggenotypes

#ENV TRIMMOMATIC /app/Trimmomatic-0.38/trimmomatic-0.38.jar



WORKDIR /out
CMD ["bash"]
