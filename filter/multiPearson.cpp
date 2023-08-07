#include "multiPearson.hpp"

typeIV::typeIV(double mean, double std, double skew, double kurt)
{
    double b1, b2, sqrt_b1, c;

    b1 = skew * skew;
    b2 = kurt;

    sqrt_b1 = skew;

    r = 6 * (b2 - b1 - 1) / (2 * b2 - 3 * b1 - 6);

    c = sqrt(16 * (r - 1) - b1 * (r - 2) * (r - 2));

    a = std * c / 4;

    gamma = -r * (r - 2) * sqrt_b1 / c;

    if (skew > 0)
    {
        a = -a;
        gamma = -gamma;
    }

    lambda = mean - (r - 2) * sqrt_b1 * std / 4;
}

multiStatePearsonator::multiStatePearsonator(const VectorXd &means,
                                             const VectorXd &stds,
                                             const VectorXd &skews,
                                             const VectorXd &kurts)
{
    for (size_t i = 0; i < means.size(); ++i)
    {
        typeIV state(means(i), stds(i), skews(i), kurts(i));
        states.push_back(state);
    }
}
