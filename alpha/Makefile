# cvector.o:
# 	g++ cvector.hpp -c cvector

all:
	g++ -std=c++17 -pthread -lfftw3 -O3 experiment.cpp -o experiment	
# -fsanitize=thread

stt:
	g++ -std=c++17 -pthread -lfftw3 -O3 stt.cpp -o stt-exp

python-ubuntu:
	g++ -O3 -Wall -shared -std=c++17 -pthread -lfftw3 -fPIC `python3 -m pybind11 --includes` cmtj.cpp -o cmtj`python3-config --extension-suffix`

python-macos:
	g++ -O3 -Wall -shared -std=c++17 -pthread -lfftw3  -undefined dynamic_lookup `python3 -m pybind11 --includes` cmtj.cpp -o cmtj`python3-config --extension-suffix`