#include <gtest/gtest.h>
#include "Matrix.h"

// Diag

TEST(Diag, Minimal) {
  Matrix<double> a(2,2,1);

  Matrix<double> result("1,0;0,1");

  EXPECT_EQ(a.diag(), result);
}

TEST(Diag, Large) {
    Matrix<double> a(100,100,1);
    a = a.diag();

    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            ASSERT_EQ(a(i,j), i == j ? 1 : 0);
        }
    }
}

TEST(Diag, OtherDiagonals) {
    Matrix<double> a(3,3,1);

    EXPECT_EQ(a.diag(1),  Matrix<double>("0,1,0;0,0,1;0,0,0"));
    EXPECT_EQ(a.diag(-1), Matrix<double>("0,0,0;1,0,0;0,1,0"));
    EXPECT_EQ(a.diag(3), Matrix<double>(3,3,0));
    EXPECT_EQ(a.diag(-3), Matrix<double>(3,3,0));
}

// Tril && Triu

TEST(GetTriangle, Triu) {
    Matrix<double> a(3,3,1);

    EXPECT_EQ(a.triu(),  Matrix<double>("1,1,1;0,1,1;0,0,1"));
    EXPECT_EQ(a.triu(-1), Matrix<double>("1,1,1;1,1,1;0,1,1"));
    EXPECT_EQ(a.triu(1), Matrix<double>("0,1,1;0,0,1;0,0,0"));
    EXPECT_EQ(a.triu(-3), Matrix<double>(3,3,1));
    EXPECT_EQ(a.triu(3), Matrix<double>(3,3,0));
}

TEST(GetTriangle, Tril) {
    Matrix<double> a(3,3,1);

    EXPECT_EQ(a.tril(),  Matrix<double>("1,0,0;1,1,0;1,1,1"));
    EXPECT_EQ(a.tril(-1), Matrix<double>("0,0,0;1,0,0;1,1,0"));
    EXPECT_EQ(a.tril(1), Matrix<double>("1,1,0;1,1,1;1,1,1"));
    EXPECT_EQ(a.tril(-3), Matrix<double>(3,3,0));
    EXPECT_EQ(a.tril(3), Matrix<double>(3,3,1));
}

// Inv

TEST(Inverse, Diagonal) {
    Matrix<double> a(3,3,5);
    a = a.diag();
    Matrix<double> a_inv = a.inv();

    EXPECT_EQ(a, a_inv.inv());
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ASSERT_EQ(a_inv(i,j), i == j ? 1 / ((double) a(i,j)) : 0);
        }
    }
}

// Solving matrix equations

TEST(SolveEquation, DirectUnique) {
    Matrix<double> A("2,1;3,-2");
    Matrix<double> b("5;7");
    Matrix<double> result(2,1);
    result(0,0) = (double)(17/7.0);
    result(1,0) = (double)(1/7.0);

    Matrix<double> x;
    ASSERT_NO_THROW(x = A | b);
    EXPECT_EQ(x, result);
}

TEST(SolveEquation, DirectNotUnique) {
    Matrix<double> A("1,3,1;1,1,-1;3,11,5");
    Matrix<double> b("9;1;35");

    Matrix<double> x;
    EXPECT_THROW(x = A | b, MatrixEquationNoUniqueSolutionException);
    EXPECT_EQ(x, Matrix<double>());
}

TEST(SolveEquation, DirectRight) {
    Matrix<double> A("1,1,1,1;1,2,3,4;1,3,6,10;1,4,10,20");
    Matrix<double> b("4,3,2,1");
    Matrix<double> result("5,-1,0,0");

    Matrix<double> x;
    ASSERT_NO_THROW(x = b / A);
    EXPECT_EQ(x, result);
}


TEST(SolveEquation, UpperTriangular) {
    Matrix<double> A(4,4,1);
    A = A.triu();
    Matrix<double> b("4;3;2;1");
    Matrix<double> result("1;1;1;1");

    Matrix<double> x;
    ASSERT_NO_THROW(x = A | b);
    EXPECT_EQ(x, result);
}

TEST(SolveEquation, LowerTriangular) {
    Matrix<double> A(4,4,1);
    A = A.tril();
    Matrix<double> b("1;2;3;4");
    Matrix<double> result("1;1;1;1");

    Matrix<double> x;
    ASSERT_NO_THROW(x = A | b);
    EXPECT_EQ(x, result);
}