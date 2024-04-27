#include <gtest/gtest.h>
#include "Matrix.h"

// Addition

TEST(Addition, Minimal) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,1,2);
  Matrix<double> c(1,1,3);

  EXPECT_EQ(a+b, c);
  EXPECT_EQ(a+a, b);  
  EXPECT_NE(a+c, b);
}

TEST(Addition, Large) {
  Matrix<double> a(100,100,1);
  Matrix<double> b(100,100,2);
  Matrix<double> c(100,100,3);

  EXPECT_EQ(a+b, c);
  EXPECT_EQ(a+a, b);  
  EXPECT_NE(a+c, b);
}

TEST(Addition, SizeMismatchCheck) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,2,1);
  Matrix<double> c(2,1,1);

  EXPECT_THROW(a+b, MatrixSizeDisparityException);
  EXPECT_THROW(a+c, MatrixSizeDisparityException);
  EXPECT_THROW(b+c, MatrixSizeDisparityException);
}

// Substraction

TEST(Substraction, Minimal) {
  Matrix<double> a(1,1,3);
  Matrix<double> b(1,1,2);
  Matrix<double> c(1,1,1);

  EXPECT_EQ(a-b, c);
  EXPECT_EQ(b-c, c);  
  EXPECT_NE(a-a, b);
}

TEST(Substraction, Large) {
  Matrix<double> a(100,100,3);
  Matrix<double> b(100,100,2);
  Matrix<double> c(100,100,1);

  EXPECT_EQ(a-b, c);
  EXPECT_EQ(b-c, c);  
  EXPECT_NE(a-a, b);
}


TEST(Substraction, SizeMismatchCheck) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,2,1);
  Matrix<double> c(2,1,1);

  EXPECT_THROW(a-b, MatrixSizeDisparityException);
  EXPECT_THROW(a-c, MatrixSizeDisparityException);
  EXPECT_THROW(b-c, MatrixSizeDisparityException);
}

// Multiply by value

TEST(MultiplyByValue, Minimal) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,1,2);
  double c = 2, d = 4;


  EXPECT_EQ(a*c, b);
  EXPECT_EQ(b*c, a*d);  
}

TEST(MultiplyByValue, Large) {
  Matrix<double> a(100,100,1);
  Matrix<double> b(100,100,2);
  double c = 2, d = 4;


  EXPECT_EQ(a*c, b);
  EXPECT_EQ(b*c, a*d);  
}

// Divide by value

TEST(DivideByValue, Minimal) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,1,2);
  Matrix<double> c(1,1,4);
  double d = 2, e = 4;


  EXPECT_EQ(b/d, a);
  EXPECT_EQ(c/e, b/d);  
}

TEST(DivideByValue, Large) {
  Matrix<double> a(100,100,1);
  Matrix<double> b(100,100,2);
  Matrix<double> c(100,100,4);
  double d = 2, e = 4;


  EXPECT_EQ(b/d, a);
  EXPECT_EQ(c/e, b/d);  
}

// Negation

TEST(Negation, Minimal) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,1,-1);

  EXPECT_EQ(a, -b);
  EXPECT_EQ(b, -a);  
}

TEST(Negation, Large) {
  Matrix<double> a(100,100,1);
  Matrix<double> b(100,100,-1);

  EXPECT_EQ(a, -b);
  EXPECT_EQ(b, -a);  
}

// Multiplication

TEST(Multiplication, Minimal) {
  Matrix<double> a(1,1,3);
  Matrix<double> b(1,1,4);

  Matrix<double> ab = a*b;
  Matrix<double> ba = b*a;
  Matrix<double> result(1,1,12);

  EXPECT_EQ(ab, result);
  EXPECT_EQ(ba, result);  
}

TEST(Multiplication, Square) {
  Matrix<double> a(5,5,3);
  Matrix<double> b(5,5,2);

  Matrix<double> ab = a*b;
  Matrix<double> ba = b*a;
  Matrix<double> result(5,5,30);

  EXPECT_EQ(ab, result);
  EXPECT_EQ(ba, result);  
}

TEST(Multiplication, SizeDifference) {
  Matrix<double> a(10,2,2);
  Matrix<double> b(2,5,3);

  Matrix<double> ab = a*b;
  Matrix<double> result(10,5,12);
  Dimensions size ={10,5};

  ASSERT_EQ(ab.getSize(), size);
  EXPECT_EQ(ab, result);
}

TEST(Multiplication, Large) {
  Matrix<double> a(100,50,1);
  Matrix<double> b(50,100,1);

  Matrix<double> ab = a*b;
  Matrix<double> result(100,100,50);
  Dimensions size ={100,100};

  ASSERT_EQ(ab.getSize(), size);
  EXPECT_EQ(ab, result);
}

TEST(Multiplication, SizeMismatchCheck) {
  Matrix<double> a(1,1,1);
  Matrix<double> b(1,2,1);
  Matrix<double> c(2,1,1);
  Matrix<double> d(2,2,1);

  EXPECT_THROW(a*c, MatrixSizeDisparityException);
  EXPECT_THROW(a*d, MatrixSizeDisparityException);
  EXPECT_THROW(b*a, MatrixSizeDisparityException);
  EXPECT_THROW(b*b, MatrixSizeDisparityException);
  EXPECT_THROW(c*c, MatrixSizeDisparityException);
  EXPECT_THROW(c*d, MatrixSizeDisparityException);
  EXPECT_THROW(d*a, MatrixSizeDisparityException);
  EXPECT_THROW(d*b, MatrixSizeDisparityException);

  EXPECT_NO_THROW(a*a);
  EXPECT_NO_THROW(a*b);
  EXPECT_NO_THROW(b*c);
  EXPECT_NO_THROW(b*d);
  EXPECT_NO_THROW(c*a);
  EXPECT_NO_THROW(c*b);
  EXPECT_NO_THROW(d*c);
  EXPECT_NO_THROW(d*d);
}

// Transposition

TEST(Transposition, Minimal) {
  Matrix<double> a(1,2,1);
  Dimensions size ={2,1};
  Matrix<double> result(size,1);

  Matrix<double> at = a.transpose();

  ASSERT_EQ(at.getSize(), size);
  EXPECT_EQ(at, result);
}

TEST(Transposition, Example) {
  Matrix<double> a("1,2,3;4,5,6;7,8,9");
  Dimensions size ={3,3};
  Matrix<double> result("1,4,7;2,5,8;3,6,9");
  
  Matrix<double> at = a.transpose();

  ASSERT_EQ(at.getSize(), size);
  EXPECT_EQ(at, result);
}

