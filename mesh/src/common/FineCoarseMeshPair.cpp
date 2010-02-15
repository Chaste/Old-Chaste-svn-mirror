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


#include "FineCoarseMeshPair.hpp"


template<unsigned DIM>
FineCoarseMeshPair<DIM>::FineCoarseMeshPair(TetrahedralMesh<DIM,DIM>& rFineMesh, QuadraticMesh<DIM>& rCoarseMesh)
    : mrFineMesh(rFineMesh),
      mrCoarseMesh(rCoarseMesh)
{
    // compute min and max values for the fine mesh nodes
    for(unsigned j=0; j<DIM; j++)
    {
        double min = 1e200;
        double max = -1e200;
     
        for(unsigned i=0; i<mrFineMesh.GetNumNodes(); i++)
        {
            if( mrFineMesh.GetNode(i)->rGetLocation()[j] < min)
            {
                min = mrFineMesh.GetNode(i)->rGetLocation()[j];
            }

            if( mrFineMesh.GetNode(i)->rGetLocation()[j] > max)
            {
                max = mrFineMesh.GetNode(i)->rGetLocation()[j];
            }
        }
        
        mMinValuesFine(j) = min;
        mMaxValuesFine(j) = max;
    }
    
    mpFineMeshBoxCollection = NULL;        
    mCounters.resize(3,0);
    
//////// only implement if going to be needed        
//////        mIdenticalMeshes = false;
//////        if(mrFineMesh.GetNumElements()==mrCoarseMesh.GetNumElements())
//////        {
//////            if(mrFineMesh.GetNumNodes()==mrCoarseMesh.GetNumVertices())
//////            {
//////                mIdenticalMeshes = true;
//////                for(unsigned i=0; i<mrFineMesh.GetNumNodes(); i++)
//////                {
//////                    for(unsigned j=0; j<DIM; j++)
//////                    {
//////                        if(fabs(mrFineMesh.GetNode(i)->rGetLocation()[j] - mrCoarseMesh.GetNode(i)->rGetLocation()[j]) > 1e-4)
//////                        {
//////                            mIdenticalMeshes = false;
//////                            break;
//////                        }
//////                    }
//////                }
//////            }
//////        }
}


template<unsigned DIM>
FineCoarseMeshPair<DIM>::~FineCoarseMeshPair()
{
    DeleteBoxCollection();
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::SetUpBoxesOnFineMesh(double boxWidth)
{
    // set up the boxes. Use a domain which is a touch larger than the fine mesh
    c_vector<double,2*DIM> min_and_max;
    for(unsigned i=0; i<DIM; i++)
    {
        min_and_max(2*i) = mMinValuesFine(i) - 0.05*fabs(mMinValuesFine(i));
        min_and_max(2*i+1) = mMaxValuesFine(i) + 0.05*fabs(mMaxValuesFine(i));
    }
    
    mpFineMeshBoxCollection = new BoxCollection<DIM>(boxWidth, min_and_max);

    // for each element, if ANY of its nodes are physically in a box, put that element 
    // in that box        
    for(unsigned i=0; i<mrFineMesh.GetNumElements(); i++)
    {
        Element<DIM,DIM>* p_element = mrFineMesh.GetElement(i);
         
        std::set<unsigned> box_indices_each_node_this_elem;
        for(unsigned j=0; j<DIM+1; j++) // num vertices per element
        {
            Node<DIM>* p_node = p_element->GetNode(j);
            unsigned box_index = mpFineMeshBoxCollection->CalculateContainingBox(p_node);
            box_indices_each_node_this_elem.insert(box_index);
        }
        
        for(std::set<unsigned>::iterator iter = box_indices_each_node_this_elem.begin();
            iter != box_indices_each_node_this_elem.end();
            ++iter)
        {
            mpFineMeshBoxCollection->rGetBox( *iter ).AddElement(p_element);
        }
    }
}


template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeFineElementsAndWeightsForCoarseQuadPoints(GaussianQuadratureRule<DIM>& rQuadRule)
{
    if(mpFineMeshBoxCollection==NULL)
    {
        EXCEPTION("Call SetUpBoxesOnFineMesh() before ComputeFineElementsAndWeightsForCoarseQuadPoints()");
    }
    
    // get the quad point (physical) positions
    QuadraturePointsGroup<DIM> quad_point_posns(mrCoarseMesh, rQuadRule);

    // resize the elements and weights vector.
    mElementsAndWeights.resize(quad_point_posns.Size());

    for(unsigned i=0; i<quad_point_posns.Size(); i++)
    {
        //std::cout << "\r " << i << " of " << quad_point_posns.Size();
        // get the box this point is in
        unsigned box_for_this_point = mpFineMeshBoxCollection->CalculateContainingBox( quad_point_posns.Get(i) );
        
        // a chaste point version of the c-vector is needed for the GetContainingElement call.
        ChastePoint<DIM> point;
        for(unsigned j=0; j<DIM; j++)
        {
            point.rGetLocation()[j]=quad_point_posns.Get(i)[j];
        }

        std::set<unsigned> test_element_indices;

        // the elements to try (initially) are those contained in the box the point is in
        // NOTE: it is possible the point to be in an element inot 'in' this box, as it is possible
        // for all element nodes to be in different boxes.
        for(typename std::set<Element<DIM,DIM>*>::iterator elem_iter = mpFineMeshBoxCollection->rGetBox(box_for_this_point).rGetElementsContained().begin();
            elem_iter != mpFineMeshBoxCollection->rGetBox(box_for_this_point).rGetElementsContained().end();
            ++elem_iter)
        {
            test_element_indices.insert((*elem_iter)->GetIndex());
        }

        unsigned elem_index;
        c_vector<double,DIM+1> weight; 

        try
        {
            // try these elements only, initially
            elem_index = mrFineMesh.GetContainingElementIndex(point, 
                                                              false, 
                                                              test_element_indices,
                                                              true /* quit if not in test_elements */);
            weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(point);
            
            mCounters[0]++;
        }
        catch(Exception& e)
        {
            // now try all the elements, but trying the elements contained in the boxes locals to this 
            // element first 
            std::set<unsigned> test_element_indices;
            
            std::set<unsigned> local_boxes = mpFineMeshBoxCollection->GetLocalBoxes(box_for_this_point);
            for(std::set<unsigned>::iterator local_box_iter = local_boxes.begin();
                local_box_iter != local_boxes.end();
                ++local_box_iter)
            {
                for(typename std::set<Element<DIM,DIM>*>::iterator elem_iter = mpFineMeshBoxCollection->rGetBox(*local_box_iter).rGetElementsContained().begin();
                    elem_iter != mpFineMeshBoxCollection->rGetBox(*local_box_iter).rGetElementsContained().end();
                    ++elem_iter)
                {
                    test_element_indices.insert((*elem_iter)->GetIndex());
                }
            }
            
            try
            {
                elem_index = mrFineMesh.GetContainingElementIndex(point, 
                                                                  false, 
                                                                  test_element_indices,
                                                                  false /* quit if not in test_elements */);
                weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(point);

                mCounters[1]++;

            }
            catch(Exception& e)
            {
                // the point is not in ANY element, store the nearest element and corresponding weights,
                // and save some information                    
                elem_index = mrFineMesh.GetNearestElementIndex(point);
                weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(point);

                mNotInMesh.push_back(i);
                mNotInMeshNearestElementWeights.push_back(weight);
                mCounters[2]++;
            }
        }

        mElementsAndWeights[i].ElementNum = elem_index;
        mElementsAndWeights[i].Weights = weight;
    }
}

template<unsigned DIM>
void FineCoarseMeshPair<DIM>::PrintStatistics()
{
    assert(mNotInMesh.size()==mCounters[2]);
    assert(mNotInMesh.size()==mNotInMeshNearestElementWeights.size());
    std::cout << "\nFineCoarseMeshPair statistics:\n";
    std::cout << "\tNum points for which containing (fine) element was found, using box containing that point, = " << mCounters[0] << "\n";
    std::cout << "\tNum points for which containing (fine) element elsewhere = " << mCounters[1] << "\n";
    std::cout << "\tNum points for which no containing element was found in fine mesh = " << mCounters[2] << "\n";
    if(mCounters[2]>0)
    {
        std::cout << "\tIndices and weights for points for which no containing element was found:\n";
        for(unsigned i=0; i<mNotInMesh.size(); i++)
        {
            std::cout << "\t\t" << mNotInMesh[i] << ", " << mNotInMeshNearestElementWeights[i] << "\n";
        }
    }
}


template<unsigned DIM>
void FineCoarseMeshPair<DIM>::DeleteBoxCollection()
{
    if(mpFineMeshBoxCollection != NULL)
    {
        delete mpFineMeshBoxCollection;
        mpFineMeshBoxCollection = NULL;
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class FineCoarseMeshPair<1>;
template class FineCoarseMeshPair<2>;
template class FineCoarseMeshPair<3>;
