#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

ublas::vector<double> generateMultivariateNormal(boost::mt19937 &generator, ublas::matrix<double> covariance);