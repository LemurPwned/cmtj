INC=-I"../core/" -I"/opt/homebrew/include/"
LIBS=-pthread -lfftw3 -L"/opt/homebrew/lib/"
CXX=clang++
stack:
	$(CXX) -std=c++17 $(INC) $(LIBS) -O3 stack_exp.cpp -o stack_exp
