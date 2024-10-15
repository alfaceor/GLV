
all:GLV_RUK4 

GLV_RUK4:GLV_RUK4.cpp
	$(CXX) -O3 -I ~/.local/include -o GLV_RUK4 GLV_RUK4.cpp

test:test.cpp
	$(CXX) -std=c++11 test.cpp -lgtest -lgtest_main -pthread -o test_run

clean:
	rm GLV_RUK4
