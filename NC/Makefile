CXX      = g++
CXXFLAGS = -std=c++17 -Wall
INCLUDE  = -Iinclude
 
SRC = src/Matrix.cpp \
      src/SLE.cpp \
      src/GaussElimination.cpp \
      src/GaussJacobi.cpp \
      src/GaussSeidel.cpp \
      src/LUCrout.cpp \
      src/LUDoolittle.cpp \
      src/LUCholesky.cpp
 
TARGET = matrixApp
 
all: $(TARGET)
 
$(TARGET): main.cpp $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDE) main.cpp $(SRC) -o $(TARGET)
 
clean:
	rm -f $(TARGET)
 
run: all
	./$(TARGET)
 