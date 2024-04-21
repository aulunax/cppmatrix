#include <iostream>
#include "Matrix.h"

int main() {

    Matrix<double> a(20000,10000, 1);
    Matrix<double> b(10000,20000, 1);
    Matrix<double> c(20000,20000);
    c = a * b;

}