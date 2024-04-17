#include <iostream>
#include "Matrix.h"

int main() {

    Matrix<int> help("4,4,7;6,3,1;0,5,1");
    help.print();

    Matrix<int> xd("1,2,3;4,5,6;7,8,9");
    xd.print();

    try {
        Matrix<int> nice = xd + help;
        nice.print();
        nice = xd;
        nice.print();
    }
    catch (MatrixSizeDisparityException e) {
        std::cout << "this sucks\n";
    }
    Matrix<int> temp = xd;
    xd = xd - xd;
    xd.print();
    xd = std::move(temp);

    Matrix<int> res = xd * help;
    res.print();

    Matrix<int> test = Matrix<int>(5,5,1) * 6 / 2;   
    test.print();

    Matrix<int> a = Matrix<int>(1000,1000,1) * Matrix<int>(1000,1000,1);
    //a.print();
    test.diag().print();
    test.tril().print();
    test.triu().print();
    test.tril(-1).print();
    test.triu(1).print();

    Matrix<int>().print();

   


}