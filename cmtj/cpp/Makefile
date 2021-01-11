INC=-I"../core/"
LIBS=-pthread -lfftw3
stt:
	g++ -std=c++17 $(INC) $(LIBS) -O3 stt.cpp -o stt-exp

vsd:
	g++ -std=c++17 $(INC) $(LIBS) -O3 vsd.cpp -o vsd

python-ubuntu:
	g++ -O3 -Wall -shared -std=c++17 -pthread -lfftw3 -fPIC `python3 -m pybind11 --includes` cmtj.cpp -o cmtj`python3-config --extension-suffix`

python-macos:
	g++ -O3 -Wall -shared -std=c++17 -pthread -lfftw3  -undefined dynamic_lookup `python3 -m pybind11 --includes` cmtj.cpp -o cmtj`python3-config --extension-suffix`