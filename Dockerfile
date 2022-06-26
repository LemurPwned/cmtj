FROM ubuntu:20.04


# let's copy all the necessary files
# we install our only 2 dependencies :) and vim for nice workflow
RUN apt-get update && \
    apt-get install -y build-essential libfftw3-dev git python3 python3-pip vim
RUN python3 -m pip install --upgrade pip
WORKDIR /scratch
# COPY cmtj /scratch/cmtj
# COPY core /scratch/core
# COPY python /scratch/python
# COPY setup.py setup.cfg ./
# RUN python3 -m pip install .
RUN git clone -b logging-fixes https://github.com/LemurPwned/cmtj.git && \
    cd cmtj && \
    python3 -m pip install numpy scipy && \
    python3 -m pip install .
WORKDIR /app
RUN python3 -c "import cmtj; from cmtj.utils import FieldScan"
