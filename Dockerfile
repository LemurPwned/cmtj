FROM ubuntu:18.04 

WORKDIR /app 

# let's copy all the necessary files
COPY alpha/*.hpp alpha/*.cpp alpha/Makefile /app/ 

# we install our only 2 dependencies :) and vim for nice workflow
RUN apt-get update && apt-get install -y build-essential libfftw3-dev \
    python3 python3-pip python-pybind11 vim && \
    python3 -m pip install pybind11

ENTRYPOINT [ "make python-ubuntu" ]