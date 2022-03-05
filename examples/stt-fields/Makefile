INC=-I"../core/" -I"/opt/homebrew/include/"
LIBS=-pthread -lfftw3 -L"/opt/homebrew/lib/"
CXX=clang++
torque:
	$(CXX) -std=c++17 $(INC) $(LIBS) -Wall -Wextra -ffast-math -O3 torque.cpp -o torque
resistance:
	$(CXX) -std=c++17 $(INC) $(LIBS) -Wall -Wextra -O3 resistance.cpp -o resistance
