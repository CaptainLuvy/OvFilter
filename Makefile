CXX = g++
CXXFLAGS = -O3
LDFLAGS = -lz -lpthread

TARGET = ovfilter
SRC = ovfilter.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)
