FROM python:3-slim-bullseye
RUN pip install biopython tqdm
COPY scripts/FastaSampler.py /usr/local/bin/FastaSampler.py
RUN chmod +x /usr/local/bin/FastaSampler.py
ENTRYPOINT ["FastaSampler.py"]

