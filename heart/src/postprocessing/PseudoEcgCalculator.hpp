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

#ifndef _PSEUDOECGALCALCULATOR_HPP_
#define _PSEUDOECGALCALCULATOR_HPP_

#include "AbstractFunctionalCalculator.hpp"
#include "TetrahedralMesh.hpp"
#include "ChastePoint.hpp"
#include "Hdf5DataReader.hpp"

/**
 * This class implements a pseudo-ECG calculator. It is a concrete class of AbstractFunctionalCalculator.
 * The pseudo-ECG is defined as the integral over the mesh of the following integrand:
 *
 * - D * grad (solution) dot grad (1/r)
 *
 *  where D is a diffusion coefficient and r is the distance between a recording electrode
 *  and a given point in the mesh.
 *
 *  References for the formula that defines the pseudo-ECG:
 *
 *  Gima K, Rudy Y Circ Res (2002) 90:889-896 (equation 1)
 *  Baher et al. Am J Physiol (2007) 292:H180-H189 (equation 5)
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class PseudoEcgCalculator : public AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    friend class TestPseudoEcgCalculator;//for testing

    Hdf5DataReader* mpDataReader; /**< An HDF5 reader from which to get the solution*/
    unsigned mNumberOfNodes; /**< Number of nodes in the mesh (got from the data reader)*/
    unsigned mNumTimeSteps;/**< Number of time steps in the simulation (got from the data reader)*/
    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& mrMesh;/**< A mesh used by the calculator*/
    ChastePoint<SPACE_DIM>& mrX; /**<The point from where we want to calculate the pseudoECG*/
    double mDiffusionCoefficient;/**<The diffusion coefficient D*/
    std::string mVariableName;/**< the variable for which we want to calculate the pseudo ecg, defaults to "V"*/

    /**
     * Get the integrand.
     *
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     */
    double GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                                c_vector<double,PROBLEM_DIM>& rU,
                                c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU);


    /**
     * Calculates the pseudo-ECG and stores the results in a vector
     *
     * @param timeStep the time step where we want to calculate the pseudoecg
     * @return the pseudo ECG at the given time step.
     *
     */
    double ComputePseudoEcgAtOneTimeStep (unsigned timeStep);

public:

    /**
     * Constructor
     *
     * @param rMesh A reference to the mesh
     * @param rX The location of the recording electrode
     * @param directory The directory where the simulation results are stored
     * @param hdf5File The file name  where the simulation results are stored
     * @param variableName  The name of the voltage variable (hardcoded to V - see todo)
     * @param makeAbsolute whether to make the path of directory absolute (using the OutputFileHandler)
     *
     *
     * ///\ todo pass in the name of the variable in the hdf5 file? it is currently hardcoded to "V".  See Doxygen above
     */
    PseudoEcgCalculator (TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                         ChastePoint<SPACE_DIM>& rX,
                         std::string directory,
                         std::string hdf5File,
                         std::string variableName = "V",
                         bool makeAbsolute = true);

    /**
     * Destructor.
     */
    ~PseudoEcgCalculator();


    /**
     * Sets the value of the diffusion coefficient (D)
     *
     * @param diffusionCoefficient The desired value of the diffusion coefficient
     *
     */
    void SetDiffusionCoefficient(double diffusionCoefficient);

    /**
     *
     * Calculates and writes the pseudo-ECG to file. the file will be named PseudoEcg.dat.
     * It will contain one column of numbers, each being the pseudoECG at each time step.
     * It will be created  by the master processor into /output relaitive to where the output
     * directory is set (by the HeartConfig or by default)
     *
     */
    void WritePseudoEcg();

};

#endif //_PSEUDOECGALCALCULATOR_HPP_


