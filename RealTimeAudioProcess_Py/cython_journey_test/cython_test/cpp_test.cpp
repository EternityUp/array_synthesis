#include "cpp_test.h"
#include <stdio.h>

Test::Test() {
  test1 = 0;
}


Test::Test(int test1) {
  this->test1 = test1;
}


Test::~Test() {}


int Test::returnFour() { return 4;}


int Test::returnFive() { return returnFour() + 1; }


Test Test::operator+(const Test& other) {
  return Test(test1 + other.test1);
}


Test Test::operator-(const Test& other) {
  return Test(test1 - other.test1);
}


MatrixXfPy Test::getMatrixXf(int d1, int d2) { 
  MatrixXfPy matrix = (MatrixXfPy)MatrixXf::Ones(d1,d2);
  matrix(0,0) = -10.0101003; // some manipulation, to show it carries over
  MatrixXfPy *ptr = new MatrixXfPy(d1, d2);
  ptr->setConstant(d1, d2, 10);
  printf("ptr->coeff(0, 0):%f\n", ptr->coeff(0, 0));
  return matrix;
}


MatrixXfPy* Test::getMatrixXfFromPtr(MatrixXfPy *ptr, int d1, int d2)
{
  ptr = new MatrixXfPy(d1, d2);
  printf("Now you are in c++ getMatrixXfFromPtr function!\n");
  ptr->setConstant(d1, d2, 10);
  printf("ptr->coeff(0, 0):%f\n", ptr->coeff(0, 0));
  return ptr;
}









