FROM ubuntu:20.04

WORKDIR /app

# let's copy all the necessary files
# we install our only 2 dependencies :) and vim for nice workflow
RUN apt-get update && \ 
    apt-get install -y build-essential libfftw3-dev git python3 python3-pip vim
RUN python3 -m pip install --upgrade pip setuptools
    
COPY cmtj /app/cmtj
COPY core /app/core 
COPY python /app/python
COPY setup.py setup.cfg ./ 
RUN python3 -m pip install .
