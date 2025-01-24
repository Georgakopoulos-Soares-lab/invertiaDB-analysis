FROM python:3.9-slim

WORKDIR /workflow

COPY mindi/coverage/utils.py /workflow
COPY requirements.txt /workflow
COPY setup.py /workflow
COPY mindi/minditool.py /workflow
COPY mindi/calc_densities.py /workflow
COPY mindi/calc_biophysical.py /workflow
COPY mindi/nonbdna_pipe.smk /workflow
COPY mindi/nonbdna_sub.sh /workflow
COPY mindi/create_design.py /workflow
COPY mindi/scheduling.py /workflow
COPY mindi/config /workflow

RUN pip install .
