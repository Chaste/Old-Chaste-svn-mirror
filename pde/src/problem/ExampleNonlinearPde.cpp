#include "ExampleNonlinearPde.hpp"

double ExampleNonlinearPde::ComputeLinearSourceTerm(ChastePoint<1> x)
{
    return 1.0;
}

double ExampleNonlinearPde::ComputeNonlinearSourceTerm(ChastePoint<1> x, double u)
{
    return 0.0;
}

c_matrix<double, 1, 1> ExampleNonlinearPde::ComputeDiffusionTerm(ChastePoint<1> , double u)
{
    return identity_matrix<double>(1)*u;
}

c_matrix<double, 1, 1> ExampleNonlinearPde::ComputeDiffusionTermPrime(ChastePoint<1> , double )
{
    return identity_matrix<double>(1)*0.0;
}

double ExampleNonlinearPde::ComputeNonlinearSourceTermPrime(ChastePoint<1> x, double u)
{
    return 0.0;
}
