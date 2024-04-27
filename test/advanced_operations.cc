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