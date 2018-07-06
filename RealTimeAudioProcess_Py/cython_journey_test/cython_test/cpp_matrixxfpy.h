#pragma once

#include "/home/xzc/eigen3.3.4/Eigen/Dense"
using namespace Eigen;
 
class MatrixXfPy : public MatrixXf { 
    public: 
        MatrixXfPy() : MatrixXf() { }
        MatrixXfPy(int rows,int cols) : MatrixXf(rows,cols) { }
        MatrixXfPy(const MatrixXf& other) : MatrixXf(other) { } 
};