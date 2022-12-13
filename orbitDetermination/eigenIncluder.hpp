
#ifndef __EIGEN_INCLUDER_HPP__
#define __EIGEN_INCLUDER_HPP__

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/OrderingMethods>
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
using Vector6d = Eigen::Vector<double, 6>;
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

template <typename Type, int Size>
using Vector = Matrix<Type, Size, 1>;

#endif
