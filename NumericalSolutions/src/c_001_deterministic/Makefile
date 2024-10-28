
CXXFLAGS = -std=c++17
# -std=c++11
all:GLV_RUK4 test utils.o test_002

utils.o:utils.h utils.cpp
	$(CXX) $(CXXFLAGS) -c utils.cpp

GLV_RUK4:GLV_RUK4.cpp utils.o
	$(CXX) $(CXXFLAGS) -O3 -I ~/.local/include -o GLV_RUK4 GLV_RUK4.cpp utils.o

test:test.cpp
	$(CXX) $(CXXFLAGS) test.cpp -lgtest -lgtest_main -pthread -o test_run

test_002:test_002.cpp
	$(CXX) $(CXXFLAGS) test_002.cpp -lgtest -lgtest_main -pthread -o test_002

clean:
	rm GLV_RUK4
