# CXX      = g++
# CXXFLAGS = -std=c++17 -Wall
# INCLUDE  = -Iinclude
 
# SRC = src/Matrix.cpp \
#       src/SLE.cpp \
#       src/GaussElimination.cpp \
#       src/GaussJacobi.cpp \
#       src/GaussSeidel.cpp \
#       src/LUCrout.cpp \
#       src/LUDoolittle.cpp \
#       src/LUCholesky.cpp
 
# TARGET = matrixApp
 
# all: $(TARGET)
 
# $(TARGET): main.cpp $(SRC)
# 	$(CXX) $(CXXFLAGS) $(INCLUDE) main.cpp $(SRC) -o $(TARGET)
 
# clean:
# 	rm -f $(TARGET)
 
# run: all
# 	./$(TARGET)
 
 CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

TARGET   = matrixApp

SRCS     = main.cpp \
           src/Matrix.cpp \
           src/SLE.cpp \
           src/GaussElimination.cpp \
           src/GaussJacobi.cpp \
           src/GaussSeidel.cpp \
           src/LUCrout.cpp \
           src/LUDoolittle.cpp \
           src/LUCholesky.cpp \
           src/EigenValue.cpp

OBJS     = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean