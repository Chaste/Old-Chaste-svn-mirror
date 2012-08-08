/*

Copyright (C) University of Oxford, 2005-2010

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

#include "LinearDiscriminantAnalysis.hpp"

void LinearDiscriminantAnalysis::InvertMatrix(const matrix<double>& rInput, matrix<double>& rInverse)
{
    // create a working copy of the input
    matrix<double> A(rInput);

    // create a permutation matrix for the LU-factorization
    permutation_matrix<std::size_t> pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);

    if( res != 0 )
    {
        EXCEPTION("Error in InvertMatrix: Matrix is singular.");
    }

    // create identity matrix of "inverse"
    rInverse = identity_matrix<double>(A.size1());

    // backsubstitute to get the inverse
    lu_substitute(A, pm, rInverse);
}

/**
 * Get the Covariance Matrix of a given input matrix
 *
 * @param rInput  The matrix to use
 * @param rInverse  The matrix to be populated with the covariant (must be the correct size already).
 */
void LinearDiscriminantAnalysis::CalculateCovariance(const matrix<double>& rInput, matrix<double>& rCov)
{
    // create a working copy of the input
    matrix<double> A(rInput);
    unsigned m = A.size1();

    for (unsigned j=0; j<A.size2(); j++)
    {
        double sum_of_column = 0.0;
        // Calculate the sum of this column
        for (unsigned i=0; i<A.size1(); i++)
        {
            sum_of_column += rInput (i,j);
        }
        // Take away the mean of the column from each entry in that column.
        for (unsigned i=0; i<A.size1(); i++)
        {
            A(i,j) = rInput(i,j) - sum_of_column/((double)(m));
        }
    }
    rCov = prod( trans(A), A );
    rCov /= (double)(m-1); // Divide by (m-1) to give an 'unbiased' estiamtor. If you want the biased one divide by m instead.
}

void LinearDiscriminantAnalysis::CalculatePooledCovariance(const std::vector<matrix<double> >& rTraining, std::vector<matrix<double> >& rCovMats, matrix<double>& rPooledCov)
{
    unsigned K = rTraining.size(); // Number of categories
    // We wipe the last two inputs so they just need to be specified as containers of the correct type.
    rPooledCov = zero_matrix<double> (rTraining[0].size2(),rTraining[0].size2());
    rCovMats.clear();

    vector<unsigned> NK (K); // Set up a vector for the number of points in each category.
    // Fill it up
    for (unsigned i=0; i<K ; i++)
    {
        NK(i) = rTraining[i].size1(); // Number of rows in each

        // Get the covariance matrix for this training set
        matrix<double> covariance_matrix (rTraining[i].size2(),rTraining[i].size2());
        CalculateCovariance(rTraining[i],covariance_matrix);
        rCovMats.push_back(covariance_matrix); // Record it.

        // Add to the pooled covariance calculation...
        rPooledCov = rPooledCov + covariance_matrix*((double)(NK (i)-1u));
    }
    // Scale the final result.
    rPooledCov = rPooledCov/((double)(sum(NK) - K));
}

void LinearDiscriminantAnalysis::CalculateMeanPoints(const std::vector<matrix<double> >& rTraining, std::vector<vector<double> >& rMeans)
{
    rMeans.clear();
    for (unsigned i=0; i<rTraining.size(); i++)
    {
        vector<double> sum_of_group = zero_vector<double> (rTraining[0].size2());
        for (unsigned j=0; j<rTraining[i].size1(); j++)
        {
            for (unsigned k=0; k<rTraining[i].size2(); k++)
            {
                sum_of_group(k) += rTraining[i](j,k);
            }
        }
        rMeans.push_back( sum_of_group/((double)(rTraining[i].size1())) );
    }
}

LinearDiscriminantAnalysis::LinearDiscriminantAnalysis(const std::vector<matrix<double> >& rTraining, bool testing)
    : mTraining(rTraining)
{
    if (!testing)
    {
        mDimension = rTraining[0].size2();
        // Check the dataset is consistent
        for (unsigned i=0; i<mTraining.size(); i++)
        {
            if (mTraining[i].size2() != mDimension)
            {
                EXCEPTION("All of the training data points must be of the same dimension.");
            }
        }

        // Set up everything we need for the discriminant functions
        CalculateMeanPoints(mTraining, mMeanTrainingPoints);
        CalculatePooledCovariance(mTraining, mCovarianceMatrices, mPooledCovarianceMatrix);

        // Do a few common-sense checks on the sizes of things.
        assert(mMeanTrainingPoints.size()==mTraining.size());
        assert(mCovarianceMatrices.size()==mTraining.size());

        // Calculate and store Sigma^-1 * mu_K for each K (reduce work needed later).
        matrix<double> inverse_pooled;
        InvertMatrix(mPooledCovarianceMatrix,inverse_pooled);
        mInvPooledDotMean.clear();
        for (unsigned K=0; K<mTraining.size(); K++)
        {
            mInvPooledDotMean.push_back(prod(inverse_pooled, mMeanTrainingPoints[K]));
        }
    }
}

std::vector<vector<double> > LinearDiscriminantAnalysis::GetMeanTrainingPoints()
{
    return mMeanTrainingPoints;
}

matrix<double> LinearDiscriminantAnalysis::GetPooledCovarianceMatrix()
{
    return mPooledCovarianceMatrix;
}

std::vector<matrix<double> > LinearDiscriminantAnalysis::GetCovarianceMatrices()
{
    return mCovarianceMatrices;
}

unsigned LinearDiscriminantAnalysis::ClassifyThisPoint(const vector<double>& rPoint)
{
    if (rPoint.size() != mDimension)
    {
        EXCEPTION("This point is not of the same dimension as the training data.");
    }

    vector<double> discrim_scores = zero_vector<double> (mTraining.size());
    double max_entry = -DBL_MAX;
    unsigned max_index = UINT_MAX;
    for (unsigned i=0; i<mTraining.size(); i++)
    {
        discrim_scores(i) = inner_prod(trans(rPoint),mInvPooledDotMean[i])
                          - 0.5*(inner_prod(trans(mMeanTrainingPoints[i]),mInvPooledDotMean[i]))
                          + log(1.0/mTraining.size());
        if (discrim_scores(i) > max_entry)
        {
            max_entry = discrim_scores(i);
            max_index = i;
        }
    }
    assert(max_index!=UINT_MAX);
    return max_index;
}
