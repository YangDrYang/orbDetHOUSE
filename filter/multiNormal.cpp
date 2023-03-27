#include "multiNormal.hpp"

ublas::vector<double> generateMultivariateNormal(boost::mt19937 &generator, ublas::matrix<double> covariance)
{
    boost::normal_distribution<double> distribution(0.0, 1.0);
    ublas::vector<double> randomVec(covariance.size1());
    for (unsigned int i = 0; i < covariance.size1(); i++)
    {
        double sum = 0;
        for (unsigned int j = 0; j < covariance.size2(); j++)
        {
            sum += distribution(generator) * sqrt(covariance(i, j));
        }
        randomVec(i) = sum;
    }
    return randomVec;
}