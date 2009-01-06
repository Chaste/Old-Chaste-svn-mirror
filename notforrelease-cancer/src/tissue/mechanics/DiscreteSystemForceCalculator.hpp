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
#ifndef DISCRETESYSTEMFORCECALCULATOR_HPP_
#define DISCRETESYSTEMFORCECALCULATOR_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "MeshBasedTissue.hpp"
#include "OutputFileHandler.hpp"

/**
 * A class for calculating the force and stress on each node in a mesh-based tissue.
 */
class DiscreteSystemForceCalculator
{
    friend class TestDiscreteSystemForceCalculator;

private:

    /** 
     * Reference to tissue.
     */
    MeshBasedTissue<2>& mrTissue;
    
    /** The mechanics used to determine the new location of the cells */
    std::vector<AbstractTwoBodyInteractionForce<2>*> mForceCollection;
    
    /**
     * Small parameter, used in GetSamplingAngles().
     */
    double mEpsilon;

    /** The file that the results of CalculateExtremalNormalForces */
    out_stream mpStressResultsFile;

    /**
     * Given a node index, returns the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);

    /**
     * Given a node index and angle of intersecting line in the range (-pi,pi],
     * returns the tangential and normal forces.
     */
    std::vector<double> CalculateFtAndFn(unsigned index, double theta);

    /**
     * Given a node index, returns a vector of sampling angles in the range (-pi,pi]
     * that can be used by GetExtremalAngles() to find the locations of local extrema
     * of the normal force.
     */
    std::vector<double> GetSamplingAngles(unsigned index);

    /**
     * Given a node index and two sampling angles, finds the location of
     * the root of the tangential force in the interval between the two
     * angles. There is no guarantee that this will lie in (-pi,pi].
     */
    double GetLocalExtremum(unsigned index, double angle1, double angle2);

    /**
     * Given a vector of sampling angles in the range (-pi,pi], returns a vector
     * of extremal angles, i.e. angles at which local extrema of the normal force
     * occur, again in the range (-pi,pi].
     */
    std::vector<double> GetExtremalAngles(unsigned index, std::vector<double> samplingAngles);

public:

    /** 
     * Constructor.
     * 
     * @param rTissue reference to the tissue
     * @param forceCollection vector of force laws present
     */
    DiscreteSystemForceCalculator(MeshBasedTissue<2>& rTissue, std::vector<AbstractTwoBodyInteractionForce<2>*> forceCollection);

    /**
     * @return the extremal normal forces on each node in the tissue.
     */
    std::vector< std::vector<double> > CalculateExtremalNormalForces();

    /**
     * Write results to file.
     * 
     * @param simulationOutputDirectory the output directory, relative to where Chaste output is stored
     */
    void WriteResultsToFile(std::string simulationOutputDirectory);

};
#endif /*DISCRETESYSTEMFORCECALCULATOR_HPP_*/
