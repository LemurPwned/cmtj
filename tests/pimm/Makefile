INC=-I"../core/" -I"/opt/homebrew/include/"
LIBS=-pthread -lfftw3 -L"/opt/homebrew/lib/"
CXX=clang++
pimm:
	$(CXX) -std=c++17 $(INC) $(LIBS) -O3 pimm.cpp -o pimm