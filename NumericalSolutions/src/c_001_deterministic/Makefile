
all:GLV_RUK4 test utils.o

utils.o:utils.h utils.cpp
	$(CXX) -c utils.cpp

GLV_RUK4:GLV_RUK4.cpp utils.o
	$(CXX) -O3 -I ~/.local/include -o GLV_RUK4 GLV_RUK4.cpp utils.o

test:test.cpp
	$(CXX) -std=c++11 test.cpp -lgtest -lgtest_main -pthread -o test_run

clean:
	rm GLV_RUK4
