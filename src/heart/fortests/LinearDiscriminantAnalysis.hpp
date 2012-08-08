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

#ifndef LINEARDISCRIMINANTANALYSIS_HPP_
#define LINEARDISCRIMINANTANALYSIS_HPP_

#include "Exception.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

class LinearDiscriminantAnalysis
{
private:
    /** Needed for serialization. */
    friend class TestLinearDiscriminantAnalysis;

    /**
     * A standard vector of matrices of training data for each group
     *
     * rows of matrices correspond to training points,
     * columns correspond to discriminant variables.
     */
    const std::vector<matrix<double> > mTraining;

    /** The dimension of the classification */
    unsigned mDimension;

    std::vector<vector<double> > mMeanTrainingPoints;

    std::vector<matrix<double> > mCovarianceMatrices;

    matrix<double> mPooledCovarianceMatrix;

    std::vector<vector<double> > mInvPooledDotMean;

    /**
     * Matrix inversion routine.
     * Uses lu_factorize and lu_substitute in uBLAS to invert a matrix.
     *
     * This is based upon Numerical Recipes in C...
     *
     * @param rInput  The matrix to be inverted
     * @param rInverse  The matrix to be populated with the inverse.
     */
    void InvertMatrix(const matrix<double>& rInput, matrix<double>& rInverse);

    /**
     * Get the Covariance Matrix of a given input matrix
     *
     * @param rInput  The matrix to use
     * @param rInverse  The matrix to be populated with the covariant (must be the correct size already).
     */
    void CalculateCovariance(const matrix<double>& rInput, matrix<double>& rCov);

    /**
     * @param rTraining  The training data, each entry in the vector is a different group.
     * @param rCovMats  Covariance matrices for each group (empty - filled by this method)
     * @param rPooledCov  The pooled covariance matrix (empty - filled by this method).
     */
    void CalculatePooledCovariance(const std::vector<matrix<double> >& rTraining,
                          std::vector<matrix<double> >& rCovMats,
                          matrix<double>& rPooledCov);

    /**
     * Calculate the mean point in each of the training groups.
     *
     * @param rTraining  The training data (mTraining is put in here by internal methods).
     * @param rMeans  Gets filled up by this method
     */
    void CalculateMeanPoints(const std::vector<matrix<double> >& rTraining, std::vector<vector<double> >& rMeans);

public:
    /**
     * Normal constructor
     *
     * @param rTraining  a vector of training data. Each entry corresponds to one of the training groups.
     * @param testing  Whether we are conducting testing on the methods and do not wish to run constructor properly.
     */
    LinearDiscriminantAnalysis(const std::vector<matrix<double> >& rTraining, bool testing=false);

    /**
     * @return  The mean point in each of the training groups.
     */
    std::vector<vector<double> > GetMeanTrainingPoints();

    /**
     * @return  The pooled covariance matrix of the training data.
     */
    matrix<double> GetPooledCovarianceMatrix();

    /**
     * @return  The covariance matrices of each training group.
     */
    std::vector<matrix<double> > GetCovarianceMatrices();

    /**
     * Perform Linear Discriminant Analysis.
     *
     * @param rPoint  The point to classify.
     * @return The index of the training group to which rPoint has been assigned.
     */
    unsigned ClassifyThisPoint(const vector<double>& rPoint);

};

#endif /* LINEARDISCRIMINANTANALYSIS_HPP_ */
