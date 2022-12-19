
#ifndef __EIGEN_INCLUDER_HPP__
#define __EIGEN_INCLUDER_HPP__

#include "../eigen-3.4.0/Eigen/Core"
#include "../eigen-3.4.0/Eigen/Dense"
#include "../eigen-3.4.0/Eigen/Sparse"
#include "../eigen-3.4.0/Eigen/SparseCholesky"
#include "../eigen-3.4.0/Eigen/SparseQR"
#include "../eigen-3.4.0/Eigen/LU"
#include "../eigen-3.4.0/Eigen/Cholesky"
#include "../eigen-3.4.0/Eigen/Geometry"
#include "../eigen-3.4.0/Eigen/OrderingMethods"
using Eigen::COLAMDOrdering;
using Eigen::LDLT;
using Eigen::LLT;
using Eigen::Matrix;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::PartialPivLU;
using Eigen::SimplicialLDLT;
using Eigen::SimplicialLLT;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::Map;
using Eigen::MatrixXi;
using Eigen::Quaterniond;
using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::SparseVector;
using Eigen::Triplet;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::placeholders::all;
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

template <typename Type, int Size>
using Vector = Matrix<Type, Size, 1>;

#endif
