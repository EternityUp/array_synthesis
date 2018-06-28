#pragma once
#ifndef AMY_MODULES_AUDIO_PROCESSING_EIGEN_WRAPPER_H
#define AMY_MODULES_AUDIO_PROCESSING_EIGEN_WRAPPER_H

#include <Eigen/Dense>

#include <complex>

using namespace Eigen;


typedef Matrix<double, Dynamic, 1> ColVectorXd;
typedef Matrix<double, 1, Dynamic> RowVectorXd;
typedef Matrix<float, Dynamic, 1> ColVectorXf;
typedef Matrix<float, 1, Dynamic> RowVectorXf;
typedef Matrix<int, Dynamic, 1> ColVectorXi;
typedef Matrix<int, 1, Dynamic> RowVectorXi;
typedef Matrix<bool, Dynamic, 1> ColVectorXb;
typedef Matrix<bool, 1, Dynamic> RowVectorXb;
typedef Matrix<double, 3, 1> ColVector3d;
typedef Matrix<double, 1, 3> RowVector3d;


typedef Matrix<std::complex<double>, Dynamic, 1> ColVectorXcd;
typedef Matrix<std::complex<double>, 1, Dynamic> RowVectorXcd;

typedef Matrix<std::complex<float>, Dynamic, 1> ColVectorXcf;
typedef Matrix<std::complex<float>, 1, Dynamic> RowVectorXcf;


typedef Matrix<float, Dynamic, Dynamic> MatrixXf;
typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
typedef Matrix<int, Dynamic, Dynamic> MatrixXi;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;


typedef Matrix <std::complex<float>, Dynamic, Dynamic> MatrixXcf;
typedef Matrix <std::complex<double>, Dynamic, Dynamic> MatrixXcd;

#endif







