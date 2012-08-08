/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef _TESTLINEARDISCRIMINANTANALYSIS_HPP_
#define _TESTLINEARDISCRIMINANTANALYSIS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include "LinearDiscriminantAnalysis.hpp"

/**
 * This file tests the linear discriminant analysis within C++.
 * This code is used by the executable application TorsadePredict in the apps folder.
 *
 * (Figures in the paper were generated with Matlab code (provided in matlab_analysis folder))
 *
 */
class TestLinearDiscriminantAnalysis : public CxxTest::TestSuite
{
private:
    matrix<double> MakeASimpleMatrix(void) throw(Exception)
    {
        // Set up a singular matrix = [0 1 2; 3 4 5; 6 7 8];
        matrix<double> singular_matrix (3, 3);
        for (unsigned i = 0; i < singular_matrix.size1 (); ++ i)
        {
            for (unsigned j = 0; j < singular_matrix.size2 (); ++ j)
            {
                singular_matrix (i, j) = 3 * i + j;
            }
        }
        return singular_matrix;
    }

public:
    void TestInvertingAMatrix(void) throw(Exception)
    {
        matrix<double> singular_matrix = MakeASimpleMatrix();
        matrix<double> nice_matrix(singular_matrix);
        nice_matrix (0,0) = 2.0; // Make this matrix nice = [2 1 2; 3 4 5; 6 7 8];

        // Set up an inverse of the correct size...
        matrix<double> inverse (3, 3);

        std::vector<matrix<double> > training;
        LinearDiscriminantAnalysis lda(training,true);
        // Attempt to invert the singular matrix.
        TS_ASSERT_THROWS_THIS(lda.InvertMatrix(singular_matrix,inverse), "Error in InvertMatrix: Matrix is singular.");

        // Invert the nice matrix
        lda.InvertMatrix(nice_matrix,inverse);

        // There are handy print methods pre-programmed on the ublas::matrix:
        //std::cout << singular_matrix << std::endl;
        //std::cout << nice_matrix << std::endl;
        //std::cout << inverse << std::endl << std::flush;

        // Check against matlab answers...
        TS_ASSERT_DELTA(inverse (0,0), 0.5, 1e-12);
        TS_ASSERT_DELTA(inverse (0,1), -1.0, 1e-12);
        TS_ASSERT_DELTA(inverse (0,2), 0.5, 1e-12);
        TS_ASSERT_DELTA(inverse (1,0), -1.0, 1e-12);
        TS_ASSERT_DELTA(inverse (1,1), -2.0/3.0, 1e-12);
        TS_ASSERT_DELTA(inverse (1,2), 2.0/3.0, 1e-12);
        TS_ASSERT_DELTA(inverse (2,0), 0.5, 1e-12);
        TS_ASSERT_DELTA(inverse (2,1), 1.0+1.0/3.0, 1e-12);
        TS_ASSERT_DELTA(inverse (2,2), -0.8 - 0.1/3.0, 1e-12);

        // Check for a ridiculously simple matrix (1D!)
        matrix<double> oned_matrix (1,1);
        oned_matrix (0,0) = 4.0;
        matrix<double> oned_inverse(oned_matrix);
        // Invert the nice matrix
        lda.InvertMatrix(oned_matrix,oned_inverse);
        TS_ASSERT_DELTA(oned_inverse (0,0), 1.0/4.0, 1e-12);
    }

    void TestCovarianceMatrix(void) throw(Exception)
    {
        matrix<double> a_matrix = MakeASimpleMatrix();
        a_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];
        a_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];

        std::vector<matrix<double> > training;
        LinearDiscriminantAnalysis lda(training,true);

        // Set up a covariant of the correct size...
        matrix<double> covariant (3, 3);
        lda.CalculateCovariance(a_matrix,covariant);

        //std::cout << a_matrix << std::endl;
        //std::cout << covariant << std::endl;

        // Check against matlab answers...
        TS_ASSERT_DELTA(covariant (0,0), 3, 1e-12);
        TS_ASSERT_DELTA(covariant (0,1), 0, 1e-12);
        TS_ASSERT_DELTA(covariant (0,2),-3, 1e-12);
        TS_ASSERT_DELTA(covariant (1,0), 0, 1e-12);
        TS_ASSERT_DELTA(covariant (1,1), 9, 1e-12);
        TS_ASSERT_DELTA(covariant (1,2), 0, 1e-12);
        TS_ASSERT_DELTA(covariant (2,0),-3, 1e-12);
        TS_ASSERT_DELTA(covariant (2,1), 0, 1e-12);
        TS_ASSERT_DELTA(covariant (2,2), 3, 1e-12);
    }

    void TestPooledCovarianceMatrix(void) throw(Exception)
    {
        // Set up the matrix from before
        matrix<double> a_matrix = MakeASimpleMatrix();
        a_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];
        a_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];

        matrix<double> b_matrix (4, 3); // Make a slightly longer matrix.
        for (unsigned i = 0; i < b_matrix.size1 (); ++ i)
        {
            for (unsigned j = 0; j < b_matrix.size2 (); ++ j)
            {
                b_matrix (i, j) = 3 * i + j;
            }
        }
        b_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (3,1) = 11.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (3,2) = 10.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];

        std::vector<matrix<double> > training_matrices;
        training_matrices.push_back(a_matrix);
        training_matrices.push_back(b_matrix);

        // We need to allocate some memory for all of the various bits and bobs...
        std::vector<matrix<double> > covariance_matrices;
        matrix<double> pooled_covariance (b_matrix.size2(),b_matrix.size2()); // Should be 3x3.

        LinearDiscriminantAnalysis lda(training_matrices,true);
        lda.CalculatePooledCovariance(training_matrices, covariance_matrices, pooled_covariance);

        TS_ASSERT_EQUALS(covariance_matrices.size(), 2u);

        matrix<double> cov1 = covariance_matrices[0];
        matrix<double> cov2 = covariance_matrices[1];

        // Check against matlab answers...
        TS_ASSERT_DELTA(cov1 (0,0), 3, 1e-12);
        TS_ASSERT_DELTA(cov1 (0,1), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (0,2),-3, 1e-12);
        TS_ASSERT_DELTA(cov1 (1,0), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (1,1), 9, 1e-12);
        TS_ASSERT_DELTA(cov1 (1,2), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (2,0),-3, 1e-12);
        TS_ASSERT_DELTA(cov1 (2,1), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (2,2), 3, 1e-12);

        // Check against matlab answers...
        TS_ASSERT_DELTA(cov2 (0,0), 18, 1e-12);
        TS_ASSERT_DELTA(cov2 (0,1), 14, 1e-12);
        TS_ASSERT_DELTA(cov2 (0,2), 4 , 1e-12);
        TS_ASSERT_DELTA(cov2 (1,0), 14, 1e-12);
        TS_ASSERT_DELTA(cov2 (1,1), 18.25, 1e-12);
        TS_ASSERT_DELTA(cov2 (1,2), 5.25, 1e-12);
        TS_ASSERT_DELTA(cov2 (2,0), 4, 1e-12);
        TS_ASSERT_DELTA(cov2 (2,1), 5.25, 1e-12);
        TS_ASSERT_DELTA(cov2 (2,2), 4.25, 1e-12);

        // Check against matlab answers...
        TS_ASSERT_DELTA(pooled_covariance (0,0), 12.0, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (0,1),  8.4, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (0,2),  1.2, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (1,0),  8.4, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (1,1),14.55, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (1,2), 3.15, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (2,0),  1.2, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (2,1), 3.15, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (2,2), 3.75, 1e-12);
    }

    void TestCalculateMeanPoints() throw(Exception)
    {
        // Set up the matrices from before
        matrix<double> a_matrix = MakeASimpleMatrix();
        a_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];
        a_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];

        matrix<double> b_matrix (4, 3); // Make a slightly longer matrix.
        for (unsigned i = 0; i < b_matrix.size1 (); ++ i)
        {
            for (unsigned j = 0; j < b_matrix.size2 (); ++ j)
            {
                b_matrix (i, j) = 3 * i + j;
            }
        }
        b_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (3,1) = 11.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (3,2) = 10.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];

        std::vector<matrix<double> > training_matrices;
        training_matrices.push_back(a_matrix);
        training_matrices.push_back(b_matrix);

        std::vector<vector<double> > mean_points;
        LinearDiscriminantAnalysis lda(training_matrices,true);
        lda.CalculateMeanPoints(training_matrices,mean_points);

        vector<double> mean1 = mean_points[0];
        vector<double> mean2 = mean_points[1];

        TS_ASSERT_DELTA(mean1(0), 1.0 , 1e-12);
        TS_ASSERT_DELTA(mean1(1), 4.0 , 1e-12);
        TS_ASSERT_DELTA(mean1(2), 7.0 , 1e-12);
        TS_ASSERT_DELTA(mean2(0), 3.0 , 1e-12);
        TS_ASSERT_DELTA(mean2(1), 5.75 , 1e-12);
        TS_ASSERT_DELTA(mean2(2), 7.75 , 1e-12);
    }

    void TestLdaAsAClass() throw(Exception)
    {
        {
            // Test an exception in the constructor - all points must be the same dimension.
            matrix<double> a_matrix = identity_matrix<double> (3, 3);
            matrix<double> b_matrix = identity_matrix<double> (4, 4);

            std::vector<matrix<double> > training_matrices;
            training_matrices.push_back(a_matrix);
            training_matrices.push_back(b_matrix);

            TS_ASSERT_THROWS_THIS(LinearDiscriminantAnalysis lda(training_matrices),
                                  "All of the training data points must be of the same dimension.");
        }

        // Set up the matrices from before
        matrix<double> a_matrix = MakeASimpleMatrix();
        a_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];
        a_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8];

        matrix<double> b_matrix (4, 3); // Make a slightly longer matrix.
        for (unsigned i = 0; i < b_matrix.size1 (); ++ i)
        {
            for (unsigned j = 0; j < b_matrix.size2 (); ++ j)
            {
                b_matrix (i, j) = 3 * i + j;
            }
        }
        b_matrix (0,2) = 8.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (2,0) = 0.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (3,1) = 11.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];
        b_matrix (3,2) = 10.0; // Make this matrix interesting = [0 1 8; 3 4 5; 0 7 8; 9 11 10];

        std::vector<matrix<double> > training_matrices;
        training_matrices.push_back(a_matrix);
        training_matrices.push_back(b_matrix);

        LinearDiscriminantAnalysis lda(training_matrices);

        std::vector<vector<double> > mean_points = lda.GetMeanTrainingPoints();
        std::vector<matrix<double> > covariance_matrices = lda.GetCovarianceMatrices();
        matrix<double> pooled_covariance = lda.GetPooledCovarianceMatrix();

        vector<double> mean1 = mean_points[0];
        vector<double> mean2 = mean_points[1];

        TS_ASSERT_DELTA(mean1(0), 1.0 , 1e-12);
        TS_ASSERT_DELTA(mean1(1), 4.0 , 1e-12);
        TS_ASSERT_DELTA(mean1(2), 7.0 , 1e-12);
        TS_ASSERT_DELTA(mean2(0), 3.0 , 1e-12);
        TS_ASSERT_DELTA(mean2(1), 5.75 , 1e-12);
        TS_ASSERT_DELTA(mean2(2), 7.75 , 1e-12);

        TS_ASSERT_EQUALS(covariance_matrices.size(), 2u);

        matrix<double> cov1 = covariance_matrices[0];
        matrix<double> cov2 = covariance_matrices[1];

        // Check against matlab answers...
        TS_ASSERT_DELTA(cov1 (0,0), 3, 1e-12);
        TS_ASSERT_DELTA(cov1 (0,1), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (0,2),-3, 1e-12);
        TS_ASSERT_DELTA(cov1 (1,0), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (1,1), 9, 1e-12);
        TS_ASSERT_DELTA(cov1 (1,2), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (2,0),-3, 1e-12);
        TS_ASSERT_DELTA(cov1 (2,1), 0, 1e-12);
        TS_ASSERT_DELTA(cov1 (2,2), 3, 1e-12);

        // Check against matlab answers...
        TS_ASSERT_DELTA(cov2 (0,0), 18, 1e-12);
        TS_ASSERT_DELTA(cov2 (0,1), 14, 1e-12);
        TS_ASSERT_DELTA(cov2 (0,2), 4 , 1e-12);
        TS_ASSERT_DELTA(cov2 (1,0), 14, 1e-12);
        TS_ASSERT_DELTA(cov2 (1,1), 18.25, 1e-12);
        TS_ASSERT_DELTA(cov2 (1,2), 5.25, 1e-12);
        TS_ASSERT_DELTA(cov2 (2,0), 4, 1e-12);
        TS_ASSERT_DELTA(cov2 (2,1), 5.25, 1e-12);
        TS_ASSERT_DELTA(cov2 (2,2), 4.25, 1e-12);

        // Check against matlab answers...
        TS_ASSERT_DELTA(pooled_covariance (0,0), 12.0, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (0,1),  8.4, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (0,2),  1.2, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (1,0),  8.4, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (1,1),14.55, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (1,2), 3.15, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (2,0),  1.2, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (2,1), 3.15, 1e-12);
        TS_ASSERT_DELTA(pooled_covariance (2,2), 3.75, 1e-12);

        vector<double> dodgy_point = zero_vector<double>(2);
        TS_ASSERT_THROWS_THIS(lda.ClassifyThisPoint(dodgy_point),
                              "This point is not of the same dimension as the training data.");

        vector<double> point = zero_vector<double>(3);

        unsigned group = lda.ClassifyThisPoint(point);
        TS_ASSERT_EQUALS(group, 0u);
    }

    void TestAProperProblem() throw(Exception)
    {
        // This is as generated by Gary's how_LDA_works.m file.

        matrix<double> a_matrix = zero_matrix<double>(5,1);
        a_matrix(0,0) = -7.861190229966832e+00;
        a_matrix(1,0) = -5.612428547792683e+00;
        a_matrix(2,0) = -5.438421658443937e+00;
        a_matrix(3,0) = -6.460175596193261e+00;
        a_matrix(4,0) = -8.593750554583927e+00;

        matrix<double> b_matrix = zero_matrix<double>(8,1);
        b_matrix(0,0) = -3.127856838023729e+00;
        b_matrix(1,0) = -4.056936972600869e+00;
        b_matrix(2,0) = -1.399133693024877e+00;
        b_matrix(3,0) = -2.288078045658942e+00;
        b_matrix(4,0) = -1.593784077187014e+00;
        b_matrix(5,0) = -3.209799951726329e+00;
        b_matrix(6,0) = -2.867324355993968e-01;
        b_matrix(7,0) = -3.196933009840541e-01;

        matrix<double> c_matrix = zero_matrix<double>(6,1);
        c_matrix(0,0) = 2.128499949162260e+00;
        c_matrix(1,0) = -8.318040651819256e-01;
        c_matrix(2,0) = 4.572823591982615e-01;
        c_matrix(3,0) = 2.542507979192909e+00;
        c_matrix(4,0) = 2.321348168622806e+00;
        c_matrix(5,0) = 5.169669566189825e-01;

        matrix<double> d_matrix = zero_matrix<double>(4,1);
        d_matrix(0,0) =  6.611901830137974e+00;
        d_matrix(1,0) =  5.181884308418971e+00;
        d_matrix(2,0) =  3.436290354230350e+00;
        d_matrix(3,0) =  5.430194341073827e+00;

        std::vector<matrix<double> > training;
        training.push_back(a_matrix);
        training.push_back(b_matrix);
        training.push_back(c_matrix);
        training.push_back(d_matrix);

        LinearDiscriminantAnalysis lda(training);

        // These are the boundaries between the groups
        // as specified by the matlab classify command
        c_vector<double,3> boundaries;
        boundaries[0] = -4.4142;
        boundaries[1] = -0.423;
        boundaries[2] = 3.177;

        for (double i=-10; i<=10; i+=0.01)
        {
            vector<double> test_point(1);
            test_point(0) = i;
            unsigned category = lda.ClassifyThisPoint(test_point);
            if (i<boundaries[0])
            {
                TS_ASSERT_EQUALS(category, 0u);
            }
            else if (i>=boundaries[0] && i<boundaries[1])
            {
                TS_ASSERT_EQUALS(category, 1u);
            }
            else if (i>=boundaries[1] && i<boundaries[2])
            {
                TS_ASSERT_EQUALS(category, 2u);
            }
            else
            {
                TS_ASSERT_EQUALS(category, 3u);
            }
        }

    }

};


#endif //_TESTLINEARDISCRIMINANTANALYSIS_HPP_
