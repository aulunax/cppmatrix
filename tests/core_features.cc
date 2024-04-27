#include <gtest/gtest.h>
#include "Matrix.h"

// Constructors

TEST(Constructor, Uninitialized) {
    Matrix<double> a;
    Dimensions size = {0,0};

    ASSERT_NO_THROW(a.getSize());
    EXPECT_EQ(a.getSize(), size);
}

TEST(Constructor, ToValueMinimal) {
    Matrix<double> a(1,1,1);
    Dimensions size = {1,1};

    EXPECT_EQ(a.getSize(), size);
    EXPECT_EQ(a(0,0), 1);
}

TEST(Constructor, ToValueLarge) {
    Matrix<double> a(100,100,1);
    Dimensions size = {100,100};

    EXPECT_EQ(a.getSize(), size);
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            ASSERT_EQ(a(i, j), 1);
        }
    }
}

TEST(Constructor, FromString) {
    Matrix<double> a("1,2,3;4,5,6;7,8,9");
    Dimensions size = {3,3};

    EXPECT_EQ(a.getSize(), size);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ASSERT_EQ(a(i, j), i*3+j+1);
        }
    }
}

TEST(Constructor, Move) {
    Matrix<double> a(1,1,1);
    Matrix<double> b(std::move(a));
    Dimensions size = {1,1};
    Dimensions uninit = {0,0};

    EXPECT_EQ(a.getSize(), uninit);
    EXPECT_EQ(b.getSize(), size);
    EXPECT_EQ(b(0,0), 1);

    Matrix<double> c(std::move(b));

    EXPECT_EQ(b.getSize(), uninit);
    EXPECT_EQ(c.getSize(), size);
    EXPECT_EQ(c(0,0), 1);
}

TEST(Constructor, Copy) {
    Matrix<double> a(1,1,2);
    Matrix<double> b(a);
    Matrix<double> c(b);

    EXPECT_EQ(a,b);
    EXPECT_EQ(c,b);
    EXPECT_EQ(a,c);
}

// Equality

TEST(Equality, Uninitialized) {
    Matrix<double> a;

    EXPECT_EQ(a, Matrix<double>());
    EXPECT_NE(a, Matrix<double>(1,1));
    EXPECT_NE(a, Matrix<double>(1,1,1));
}

TEST(Equality, Minimal) {
    Matrix<double> a(1,1,1);

    EXPECT_EQ(a, Matrix<double>(1,1,1));
    EXPECT_NE(a, Matrix<double>(1,1,2));
}

TEST(Equality, Large) {
    Matrix<double> a(100,100,1);

    EXPECT_EQ(a, Matrix<double>(100,100,1));
    EXPECT_NE(a, Matrix<double>(100,100,2));
}

TEST(Equality, SizeMismatchCheck) {
    Matrix<double> a(1,1,1);

    EXPECT_NE(a, Matrix<double>(2,1,1));
    EXPECT_NE(a, Matrix<double>(1,2,1));
    EXPECT_NE(a, Matrix<double>(2,2,1));
}

// Assignment

TEST(Assignment, Move) {
    Matrix<double> a(1,1,1);
    Matrix<double> b(1,1,3);
    Matrix<double> c;
    Dimensions size = {1,1};
    Dimensions uninit = {0,0};

    b = std::move(a);
    EXPECT_EQ(a.getSize(), uninit);
    EXPECT_EQ(b.getSize(), size);
    EXPECT_EQ(b(0,0), 1);
    c = std::move(b);
    EXPECT_EQ(b.getSize(), uninit);
    EXPECT_EQ(c.getSize(), size);
    EXPECT_EQ(c(0,0), 1);
    a = std::move(c);
    EXPECT_EQ(c.getSize(), uninit);
    EXPECT_EQ(a.getSize(), size);
    EXPECT_EQ(a(0,0), 1);
}

TEST(Assignment, Copy) {
    Matrix<double> a(1,1,1);
    Matrix<double> b(1,1,3);
    Matrix<double> c;
    b = a;
    c = b;
    EXPECT_EQ(a,b);
    EXPECT_EQ(c,b);
    EXPECT_EQ(a,c);
}