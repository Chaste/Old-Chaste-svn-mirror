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

#include "HeartGeometryInformation.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

// Area of the septum considered to belong to the each ventricle (relative to 1)
template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM_SIZE = 2.0/3.0;

template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM_SIZE = 1.0/3.0;

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                              const std::string& rEpiFile,
                                                              const std::string& rEndoFile,
                                                              bool indexFromZero)
   : mpMesh(&rMesh)
{
    DistanceMapCalculator<SPACE_DIM, SPACE_DIM> distance_calculator(*mpMesh);

    // Get nodes defining each surface
    GetNodesAtSurface(rEpiFile, mEpiSurface, indexFromZero);
    GetNodesAtSurface(rEndoFile, mEndoSurface, indexFromZero);

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(mEpiSurface, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(mEndoSurface, mDistMapEndocardium);
    mNumberOfSurfacesProvided = 2;
}

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation (AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                               const std::string& rEpiFile,
                                                               const std::string& rLVFile,
                                                               const std::string& rRVFile,
                                                               bool indexFromZero)
    : mpMesh(&rMesh)
{
    DistanceMapCalculator<SPACE_DIM, SPACE_DIM> distance_calculator(*mpMesh);

    // Get nodes defining each surface

    GetNodesAtSurface(rEpiFile, mEpiSurface, indexFromZero);
    GetNodesAtSurface(rLVFile, mLVSurface, indexFromZero);
    GetNodesAtSurface(rRVFile, mRVSurface, indexFromZero);

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(mEpiSurface, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(mLVSurface, mDistMapLeftVentricle);
    distance_calculator.ComputeDistanceMap(mRVSurface, mDistMapRightVentricle);

    mNumberOfSurfacesProvided = 3;
}


template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation (std::string nodeHeterogeneityFileName)
{
    mpMesh = NULL;
    std::ifstream heterogeneity_file;
    heterogeneity_file.open(nodeHeterogeneityFileName.c_str());

    if (!(heterogeneity_file.is_open()))
    {
        heterogeneity_file.close();
        std::stringstream ss;
        ss << "Could not open heterogeneities file (" << nodeHeterogeneityFileName << ")";
        EXCEPTION(ss.str());
    }

    while(!heterogeneity_file.eof())
    {
        int value;

        heterogeneity_file >> value;

        // format error (for example read a double), or value not equal to 0, 1, or 2.
        if( (heterogeneity_file.fail() && !heterogeneity_file.eof()) || value < 0 || value > 2)
        {
            heterogeneity_file.close();
            std::stringstream ss;
            ss << "A value in the heterogeneities file (" << nodeHeterogeneityFileName
               << ") is out of range (or not an integer). It should be epi = 0, mid = 1, endo = 2";
            EXCEPTION(ss.str());
        }

        if(!heterogeneity_file.eof())
        {
            if(value==0)
            {
                mLayerForEachNode.push_back(EPI);
            }
            else if(value==1)
            {
                mLayerForEachNode.push_back(MID);
            }
            else
            {
                assert(value==2);
                mLayerForEachNode.push_back(ENDO);
            }
        }
    }

    heterogeneity_file.close();
}

template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::ProcessLine(
        const std::string& line, std::set<unsigned>& rSurfaceNodeIndexSet, unsigned offset) const
{
    std::stringstream line_stream(line);
    while (!line_stream.eof())
    {
        unsigned item;
        line_stream >> item;
        // If offset==1 then shift the nodes, since we are assuming MEMFEM format (numbered from 1 on)
        if (item == 0 && offset != 0) 
        {
            EXCEPTION("Error when reading surface file.  It was assumed not to be indexed from zero, but zero appeared in the list.");
        }
        rSurfaceNodeIndexSet.insert(item-offset);
    }
}


template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::GetNodesAtSurface(
        const std::string& surfaceFile, std::vector<unsigned>& rSurfaceNodes, bool indexFromZero) const
{
    // Open the file defining the surface
    std::ifstream file_stream;
    unsigned offset=0;
    if (indexFromZero == false)
    {
        offset=1;
    }

    file_stream.open(surfaceFile.c_str());
    if (!file_stream.is_open())
    {
        EXCEPTION("Wrong surface definition file name " + surfaceFile);
    }

    // Temporary storage for the nodes, helps discarding repeated values
    std::set<unsigned> surface_original_node_index_set;

    // Loop over all the triangles and add node indexes to the set
    std::string line;
    getline(file_stream, line);
    do
    {
        ProcessLine(line, surface_original_node_index_set, offset);

        getline(file_stream, line);
    }
    while(!file_stream.eof());
    file_stream.close();

    // Make vector big enough
    rSurfaceNodes.reserve(surface_original_node_index_set.size());

    if (mpMesh->rGetNodePermutation().empty())
    {
        // Copy the node indexes from the set to the vector as they are
        for(std::set<unsigned>::iterator node_index_it=surface_original_node_index_set.begin();
            node_index_it != surface_original_node_index_set.end();
            node_index_it++)
        {
            rSurfaceNodes.push_back(*node_index_it);
        }
    }
    else
    {
        // Copy the original node indices from the set to the vector applying the permutation
        for(std::set<unsigned>::iterator node_index_it=surface_original_node_index_set.begin();
            node_index_it != surface_original_node_index_set.end();
            node_index_it++)
        {
            rSurfaceNodes.push_back(mpMesh->rGetNodePermutation()[*node_index_it]);
        }
    }
}



template<unsigned SPACE_DIM>
HeartRegionType HeartGeometryInformation<SPACE_DIM>::GetHeartRegion(unsigned nodeIndex) const
{

    if (mDistMapRightVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
        mDistMapRightVentricle[nodeIndex] >= mDistMapLeftVentricle[nodeIndex])
    {
        return HeartRegionCode::LEFT_VENTRICLE_WALL;
    }

    if (mDistMapLeftVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
        mDistMapLeftVentricle[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
    {
        return HeartRegionCode::RIGHT_VENTRICLE_WALL;
    }

    if (mDistMapEpicardium[nodeIndex] >= mDistMapLeftVentricle[nodeIndex] &&
        mDistMapEpicardium[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
    {
        if (mDistMapLeftVentricle[nodeIndex]
            < LEFT_SEPTUM_SIZE*(mDistMapLeftVentricle[nodeIndex] + mDistMapRightVentricle[nodeIndex]))
        {
            return HeartRegionCode::LEFT_SEPTUM;
        }
        else
        {
            return HeartRegionCode::RIGHT_SEPTUM;
        }
    }

    return HeartRegionCode::UNKNOWN;
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::GetDistanceToEndo(unsigned nodeIndex)
{
    // General case where you provide 3 surfaces: LV, RV, epicardium
    if ( mNumberOfSurfacesProvided == 3)
    {
        HeartRegionType node_region = GetHeartRegion(nodeIndex);
        switch(node_region)
        {
            case HeartRegionCode::LEFT_VENTRICLE_WALL:
            case HeartRegionCode::LEFT_VENTRICLE_SURFACE:
                return mDistMapLeftVentricle[nodeIndex];
                break;

            case HeartRegionCode::RIGHT_VENTRICLE_WALL:
            case HeartRegionCode::RIGHT_VENTRICLE_SURFACE:
                return mDistMapRightVentricle[nodeIndex];
                break;

            case HeartRegionCode::LEFT_SEPTUM:
                return mDistMapLeftVentricle[nodeIndex];
                break;

            case HeartRegionCode::RIGHT_SEPTUM:
                return mDistMapRightVentricle[nodeIndex] ;
                break;

            case HeartRegionCode::UNKNOWN:
                #define COVERAGE_IGNORE
                std::cerr << "Wrong distances node: " << nodeIndex << "\t"
                          << "Epi " << mDistMapEpicardium[nodeIndex] << "\t"
                          << "RV " << mDistMapRightVentricle[nodeIndex] << "\t"
                          << "LV " << mDistMapLeftVentricle[nodeIndex]
                          << std::endl;

                // Make wall_thickness=0 as in Martin's code
                return 0.0;
                break;
                #undef COVERAGE_IGNORE

            default:
                NEVER_REACHED;
        }
    }
    // Simplified case where you only provide epi and endo surface definitions
    else
    {
        return mDistMapEndocardium[nodeIndex];
    }

    // gcc wants to see a return statement at the end of the method.
    NEVER_REACHED;
    return 0.0;
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::GetDistanceToEpi(unsigned nodeIndex)
{
    return mDistMapEpicardium[nodeIndex];
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::CalculateRelativeWallPosition(unsigned nodeIndex)
{

    double dist_endo = GetDistanceToEndo(nodeIndex);
    double dist_epi = GetDistanceToEpi(nodeIndex);

    double relative_position;

    if ( (dist_endo + dist_epi) != 0 )
    {
       relative_position = dist_endo / (dist_endo + dist_epi);
    }
    else
    {
        /*
         *  A node contained on both epicardium and lv (or rv) surfaces has wall thickness 0/0.
         *  By setting its value to 0 we consider it contained only on the lv (or rv) surface.
         */
        relative_position = 0;
    }
    return relative_position;
}

template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::DetermineLayerForEachNode(double epiFraction, double endoFraction)
{
    if (epiFraction+endoFraction>1)
    {
        EXCEPTION("The sum of fractions of epicardial and endocardial layers must be lesser than 1");
    }

    if ((endoFraction<0) || (epiFraction<0))
    {
        EXCEPTION("A fraction of a layer must be positive");
    }

    mLayerForEachNode.resize(mpMesh->GetNumNodes());
    for(unsigned i=0; i<mpMesh->GetNumNodes(); i++)
    {
        double position = CalculateRelativeWallPosition(i);
        if (position<endoFraction)
        {
            mLayerForEachNode[i] = ENDO;
        }
        else if (position<(1-epiFraction))
        {
            mLayerForEachNode[i] = MID;
        }
        else
        {
            mLayerForEachNode[i] = EPI;
        }
    }
}


template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::WriteLayerForEachNode(std::string outputDir, std::string file)
{
    OutputFileHandler handler(outputDir,false);
    if (PetscTools::AmMaster())
    {
        out_stream p_file = handler.OpenOutputFile(file);

        assert(mLayerForEachNode.size()>0);
        for(unsigned i=0; i<mpMesh->GetNumNodes(); i++)
        {
            if(mLayerForEachNode[i]==EPI)
            {
                *p_file << "0\n";
            }
            else if(mLayerForEachNode[i]==MID)
            {
                *p_file << "1\n";
            }
            else // endo
            {
                *p_file << "2\n";
            }
        }

        p_file->close();
    }
    PetscTools::Barrier("HeartGeometryInformation::WriteLayerForEachNode"); // Make other processes wait until we're done
}


template<unsigned SPACE_DIM>
ChasteCuboid<SPACE_DIM> HeartGeometryInformation<SPACE_DIM>::CalculateBoundingBoxOfSurface(const std::vector<unsigned>& rSurfaceNodes)
{

    assert(rSurfaceNodes.size()>0);
    //Set min to DBL_MAX etc.
    c_vector<double, SPACE_DIM> my_minimum_point;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        my_minimum_point[i]=DBL_MAX; //Start with max and work down to actual
    }
    c_vector<double, SPACE_DIM> my_maximum_point=-my_minimum_point;

    //Iterate through the set of points on the surface
    for (unsigned surface_index=0; surface_index<rSurfaceNodes.size(); surface_index++)
    {
        unsigned global_index=rSurfaceNodes[surface_index];
        if (mpMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(global_index) )
        {
            c_vector<double, SPACE_DIM> position = mpMesh->GetNode(global_index)->rGetLocation();
            //Update max/min
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                if (position[i] < my_minimum_point[i])
                {
                    my_minimum_point[i] = position[i];
                }
                if (position[i] > my_maximum_point[i])
                {
                    my_maximum_point[i] = position[i];
                }
            }
        }
    }

    //Share the local data and reduce over all processes
    c_vector<double, SPACE_DIM> global_minimum_point;
    c_vector<double, SPACE_DIM> global_maximum_point;
    MPI_Allreduce(&my_minimum_point[0], &global_minimum_point[0], SPACE_DIM, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(&my_maximum_point[0], &global_maximum_point[0], SPACE_DIM, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);


    ChastePoint<SPACE_DIM> min(global_minimum_point);
    ChastePoint<SPACE_DIM> max(global_maximum_point);

    return ChasteCuboid<SPACE_DIM>(min, max);
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//template class HeartGeometryInformation<1>;
template class HeartGeometryInformation<2>;
template class HeartGeometryInformation<3>;
