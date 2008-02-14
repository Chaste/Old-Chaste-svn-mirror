#include "UblasCustomFunctions.hpp"


c_vector<double, 1> Create_c_vector(double x)
{
    c_vector<double, 1> v;
    v[0] = x;
    return v;
}

c_vector<double, 2> Create_c_vector(double x, double y)
{
    c_vector<double, 2> v;
    v[0] = x;
    v[1] = y;
    return v;
}

c_vector<double, 3> Create_c_vector(double x, double y, double z)
{
    c_vector<double, 3> v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
}
