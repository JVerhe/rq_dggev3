CXX = g++
FC = gfortran

CXXFLAGS = -O2 -std=c++17 -Iinclude -I/usr/include/eigen/
LDFLAGS = -llapack -lblas

SRC_CPP = $(wildcard src/*.cpp)
SRC_F = src/dggev3_qr.f
TARGET = main

OBJ_CPP = $(SRC_CPP:.cpp=.o)
OBJ_F = $(SRC_F:.f=.o)

$(TARGET): $(OBJ_F) $(OBJ_CPP)
	$(CXX) $(OBJ_CPP) $(OBJ_F) -o $@ $(LDFLAGS)

%.o: %.f
	$(FC) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ_CPP) $(OBJ_F) $(TARGET)
