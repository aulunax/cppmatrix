#include <iostream>
#include "Matrix.h"

int main() {



    Matrix<double> diag(3,3,3);
    diag = diag.diag();
    Matrix<double> diagInv = diag.inv();
    diagInv.print();


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

    //Matrix<int> a = Matrix<int>(300,100,1) * Matrix<int>(100,500,1);
    //a.transpose();
    Matrix<int> z = Matrix<int>(2000,1000,1) * Matrix<int>(1000,2000,1);
    //z.transpose();
    test.diag().print();
    test.diag(1).print();
    test.tril().print();
    test.triu().print();
    test.tril(-1).print();
    test.triu(1).print();

    Matrix<int>().print();

    Matrix<int> c = Matrix<int>(2,2,1).tril();
    c.print();
    c = Matrix<int>(2,2,1).triu();
    c.print();
    c = Matrix<int>(1,1,1).diag();
    c.print();
    c = Matrix<int>(2,2,1).diag();
    c.print();
    c = Matrix<int>(1, 1, 1).tril();
    c.print();
    c = Matrix<int>(1, 1, 1).triu();
    c.print();

    if (Matrix<int>(1, 1, 1).triu() == Matrix<int>(1, 1, 1).tril())
        std::cout << "works\n";

    if (Matrix<int>(2, 2, 1).triu() != Matrix<int>(2, 2, 1).tril())
        std::cout << "works\n";

    (-xd).print();

    Matrix<int> tr("1,2,3;4,5,6");
    tr.print();
    tr.transpose().print();

    try {
        Matrix<double> A("2,1;3,-2");
        Matrix<double> b("5;7");

        Matrix<double> x = A | b;
        std::cout << "Solution:\n\n";
        x.print();
    }
    catch (MatrixEquationNoUniqueSolutionException e) {
        std::cout << "No unique solutions\n\n";
    }

    try {
        Matrix<double> A("1,3,1;1,1,-1;3,11,5");
        Matrix<double> b("9;1;35");

        Matrix<double> x = A | b;
        std::cout << "Solution:\n\n";
        x.print();
    }
    catch (MatrixEquationNoUniqueSolutionException e) {
        std::cout << "No unique solutions\n\n";
    }

    try {
        Matrix<double> A("1,1,1,1;1,2,3,4;1,3,6,10;1,4,10,20");
        Matrix<double> b("4,3,2,1");

        Matrix<double> x = b / A;
        std::cout << "Solution:\n\n";
        x.print();
    }
    catch (MatrixEquationNoUniqueSolutionException e) {
        std::cout << "No unique solutions\n\n";
    }

    try {
        Matrix<double> A("1,1,1,1;1,2,3,4;1,3,6,10;1,4,10,20");
        Matrix<double> b("4,3,2,1");

        Matrix<double> x = b / A;
        std::cout << "Solution:\n\n";
        x.print();
    }
    catch (MatrixEquationNoUniqueSolutionException e) {
        std::cout << "No unique solutions\n\n";
    }
}