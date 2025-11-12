CXX = g++
FC = gfortran

CXXFLAGS = -O2 -std=c++17 -Iinclude -I/usr/include/eigen/
LDFLAGS = -llapack -lblas

TARGET = dggev3_example
CPP_SOURCES = src/main.cpp src/dggev3_wrapper.cpp
F_SOURCES = src/dggev3_qr.f src/dggev3_rq.f

CPP_OBJECTS = $(CPP_SOURCES:.cpp=.o)
F_OBJECTS = $(F_SOURCES:.f=.o)

$(TARGET): $(CPP_OBJECTS) $(F_OBJECTS)
	$(CXX) $(CPP_OBJECTS) $(F_OBJECTS) -o $@ $(LDFLAGS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

src/%.o: src/%.f
	$(FC) -c $< -o $@

clean:
	rm -f $(CPP_OBJECTS) $(F_OBJECTS) $(TARGET)
