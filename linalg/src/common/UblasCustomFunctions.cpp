/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


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

c_vector<double,3> CalculateSmallestEigenvector(c_matrix<double,3,3> &A)
{
    int info;
    c_vector<double, 3 > WR;
    c_vector<double, 3 > WI;
    c_vector<double, 4*3 > WORK;
    c_matrix<double, 3, 3> VL;
    c_matrix<double, 3, 3> VR;
    
    c_vector<double, 3> output;
    
    char N = 'N';
    char V = 'V';
    int size = 3;
    int four_times_size = 4*3;
    
    c_matrix<double, 3, 3> a_transpose;
    noalias(a_transpose) = trans(A);    
    
    dgeev_(&N,&V,&size,a_transpose.data(),&size,WR.data(),WI.data(),VL.data(),&size,VR.data(),&size,WORK.data(),&four_times_size,&info);
    assert(info==0);    
    assert(norm_2(WI) == 0.0); // We've found a complex eigenvalue... 

    
    int index_of_smallest=0;    
    double min_eigenvalue = abs(WR(0));
    
    if (abs(WR(1)) < min_eigenvalue)
    {
        index_of_smallest = 1;
        min_eigenvalue = abs(WR(1));
    }

    if (abs(WR(2)) < min_eigenvalue)
    {
        index_of_smallest = 2;
        min_eigenvalue = abs(WR(2));
    }
            
    output(0) = VR(index_of_smallest,0);
    output(1) = VR(index_of_smallest,1);
    output(2) = VR(index_of_smallest,2);
    
    return output;
    
}
