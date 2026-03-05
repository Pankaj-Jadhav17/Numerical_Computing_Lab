CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra

# Collect all .cpp source files
SRCS = src/Matrix.cpp src/GaussElimination.cpp src/LU.cpp main.cpp

all:
	$(CXX) $(CXXFLAGS) $(SRCS) -Iinclude -o matrix_app

run:
	./matrix_app

clean:
	rm -f matrix_app