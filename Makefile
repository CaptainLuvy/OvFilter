CXX = g++
CXXFLAGS = -O3 -std=c++11 -pthread -march=native

all: ovfilter_cpp

ovfilter_cpp: ovfilter.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lz

clean:
	rm -f ovfilter_cpp
