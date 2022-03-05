INC=-I"../core/" -I"/opt/homebrew/include/"
LIBS=-pthread -lfftw3 -L"/opt/homebrew/lib/"
CXX=clang++
vsd:
	$(CXX) -std=c++17 $(INC) $(LIBS) -O3 vsd.cpp -o vsd
