# CXX = g++
# CXXFLAGS = -std=c++17 -Wall

# INCLUDE = -Iinclude

# SRC = src/Matrix.cpp \
#       src/SLE.cpp \
#       src/GaussElimination.cpp \
#       src/GaussJacobi.cpp \
#       src/LUCrout.cpp \
#       src/LUDoolittle.cpp \
#       src/LUCholesky.cpp

# OBJ = $(SRC:.cpp=.o)

# TARGET = matrixApp

# all: $(TARGET)

# $(TARGET): $(OBJ) main.o
# 	$(CXX) $(CXXFLAGS) $(OBJ) main.o -o $(TARGET)

# main.o: main.cpp
# 	$(CXX) $(CXXFLAGS) $(INCLUDE) -c main.cpp

# %.o: %.cpp
# 	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# clean:
# 	rm -f $(OBJ) main.o $(TARGET)


CXX      = g++
CXXFLAGS = -std=c++17 -Wall
INCLUDE  = -Iinclude
 
SRC = src/Matrix.cpp \
      src/SLE.cpp \
      src/GaussElimination.cpp \
      src/GaussJacobi.cpp \
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
 