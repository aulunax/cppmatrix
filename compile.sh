#!/bin/bash

# compiles the src/main.cpp, for test purposes only
g++ -O3 src/Matrix.cpp src/main.cpp src/MatrixMisc.cpp -I./inc  -o test
