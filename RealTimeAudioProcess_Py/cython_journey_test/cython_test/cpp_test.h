#pragma once

#include "cpp_matrixxfpy.h"

class Test {
public:
  int test1;
  Test();
  Test(int test1);
  ~Test();
  int returnFour();
  int returnFive();
  Test operator+(const Test& other);
  Test operator-(const Test& other);
  MatrixXfPy *ptr;
  MatrixXfPy getMatrixXf(int d1, int d2); 
  MatrixXfPy* getMatrixXfFromPtr(MatrixXfPy *ptr, int d1, int d2);
};


