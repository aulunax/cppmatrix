#!/bin/bash
g++ -O3 src/Matrix.cpp src/main.cpp src/MatrixMisc.cpp -I./inc  -o test
c++ -w -Wfatal-errors -O3 -Wall -shared -std=c++17 -I ./inc -fPIC $(python3 -m pybind11 --includes) src/MatrixMisc.cpp src/Matrix.cpp src/py_matrix.cpp -o py_matrix$(python3-config --extension-suffix) 
stubgen -m py_matrix -o .