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
    mMinMaxValuesInFineMesh = mrFineMesh.GetExtremes();

    mpFineMeshBoxCollection = NULL;
    mCounters.resize(4,0);

//// Bring back identical-meshes functionality if needed
//    // check whether the two meshes are the same (ie the linear part of the quad mesh is the
//    // linear mesh, by checking whether the number of elements match and the vertex indices of
//    // each element match.
//    mIdenticalMeshes = false;
//    if(mrFineMesh.GetNumElements()==mrCoarseMesh.GetNumElements())
//    {
//        mIdenticalMeshes = true;
//        for (unsigned i=0; i<mrFineMesh.GetNumElements(); i++)
//        {
//            for (unsigned j=0; j<mrFineMesh.GetElement(i)->GetNumNodes(); j++)
//            {
//                if (mrFineMesh.GetElement(i)->GetNodeGlobalIndex(j)!=mrCoarseMesh.GetElement(i)->GetNodeGlobalIndex(j) )
//                {
//                    mIdenticalMeshes = false;
//                    break;
//                }
//            }
//        }
//    }
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
    c_vector<double,2*DIM> extended_min_and_max;
    for(unsigned i=0; i<DIM; i++)
    {
        // subtract from the minima
        extended_min_and_max(2*i) = mMinMaxValuesInFineMesh(2*i) - 0.05*fabs(mMinMaxValuesInFineMesh(2*i));
        // add to the maxima
        extended_min_and_max(2*i+1) = mMinMaxValuesInFineMesh(2*i+1) + 0.05*fabs(mMinMaxValuesInFineMesh(2*i+1));
    }

    if(boxWidth < 0)
    {
        // use default value = max( max_edge_length, w20),  where w20 is the width corresponding to 
        // 20 boxes in the x-direction 

        // BoxCollection creates an extra box so divide by 19 not 20.  Add a little bit on to ensure
        // minor numerical fluctuations don't change the answer.
        boxWidth = (extended_min_and_max(1) - extended_min_and_max(0))/19.000000001;

        // determine the maximum edge length
        double max_edge_length = -1;

        for (typename TetrahedralMesh<DIM,DIM>::EdgeIterator edge_iterator = mrFineMesh.EdgesBegin();
             edge_iterator!=mrFineMesh.EdgesEnd();
             ++edge_iterator)
        {
            c_vector<double, 3> location1 = edge_iterator.GetNodeA()->rGetLocation();
            c_vector<double, 3> location2 = edge_iterator.GetNodeB()->rGetLocation();
            double edge_length = norm_2(location1-location2);

            if(edge_length>max_edge_length)
            {
                max_edge_length = edge_length;
            }
        }
        
        if(boxWidth < max_edge_length)
        {
            boxWidth = 1.1*max_edge_length;
        }
    }
 
    mpFineMeshBoxCollection = new BoxCollection<DIM>(boxWidth, extended_min_and_max);
    mpFineMeshBoxCollection->SetupAllLocalBoxes();

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
void FineCoarseMeshPair<DIM>::ComputeFineElementsAndWeightsForCoarseQuadPoints(GaussianQuadratureRule<DIM>& rQuadRule,
                                                                               bool safeMode)
{
    if(mpFineMeshBoxCollection==NULL)
    {
        EXCEPTION("Call SetUpBoxesOnFineMesh() before ComputeFineElementsAndWeightsForCoarseQuadPoints()");
    }

    // get the quad point (physical) positions
    QuadraturePointsGroup<DIM> quad_point_posns(mrCoarseMesh, rQuadRule);

    // resize the elements and weights vector.
    mFineMeshElementsAndWeights.resize(quad_point_posns.Size());

    for(unsigned i=0; i<quad_point_posns.Size(); i++)
    {
        //std::cout << "\r " << i << " of " << quad_point_posns.Size() << std::flush;
        
        // get the box this point is in
        unsigned box_for_this_point = mpFineMeshBoxCollection->CalculateContainingBox( quad_point_posns.Get(i) );

        // a chaste point version of the c-vector is needed for the GetContainingElement call.
        ChastePoint<DIM> point(quad_point_posns.Get(i));

        ComputeFineElementAndWeightForGivenPoint(point, safeMode, box_for_this_point, i);
    }
}


template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeFineElementsAndWeightsForCoarseNodes(bool safeMode)
{
    if(mpFineMeshBoxCollection==NULL)
    {
        EXCEPTION("Call SetUpBoxesOnFineMesh() before ComputeFineElementsAndWeightsForCoarseNodes()");
    }

    // resize the elements and weights vector.
    mFineMeshElementsAndWeights.resize(mrCoarseMesh.GetNumNodes());

    for(unsigned i=0; i<mrCoarseMesh.GetNumNodes(); i++)
    {
        //std::cout << "\r " << i << " of " << mrCoarseMesh.GetNumNodes() << std::flush;
        
        Node<DIM>* p_node = mrCoarseMesh.GetNode(i);
        
        // get the box this point is in
        unsigned box_for_this_point = mpFineMeshBoxCollection->CalculateContainingBox( p_node->rGetModifiableLocation() );

        // a chaste point version of the c-vector is needed for the GetContainingElement call.
        ChastePoint<DIM> point(p_node->rGetLocation());

        ComputeFineElementAndWeightForGivenPoint(point, safeMode, box_for_this_point, i);
    }
}


template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeFineElementAndWeightForGivenPoint(ChastePoint<DIM>& rPoint, 
                                                                       bool safeMode,
                                                                       unsigned boxForThisPoint,
                                                                       unsigned index)
{                                                            
    std::set<unsigned> test_element_indices;

    // the elements to try (initially) are those contained in the box the point is in
    // NOTE: it is possible the point to be in an element inot 'in' this box, as it is possible
    // for all element nodes to be in different boxes.
    for(typename std::set<Element<DIM,DIM>*>::iterator elem_iter = mpFineMeshBoxCollection->rGetBox(boxForThisPoint).rGetElementsContained().begin();
        elem_iter != mpFineMeshBoxCollection->rGetBox(boxForThisPoint).rGetElementsContained().end();
        ++elem_iter)
    {
        test_element_indices.insert((*elem_iter)->GetIndex());
    }

    unsigned elem_index;
    c_vector<double,DIM+1> weight;

    try
    {
        //std::cout << "\n" << "# test elements initially " << test_element_indices.size() << "\n";
        // try these elements only, initially
        elem_index = mrFineMesh.GetContainingElementIndex(rPoint,
                                                          false,
                                                          test_element_indices,
                                                          true /* quit if not in test_elements */);
        weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

        mCounters[0]++;
    }
    catch(Exception& e)
    {
        // now try all the elements, but trying the elements contained in the boxes local to this
        // element first
        std::set<unsigned> test_element_indices;

        std::set<unsigned> local_boxes = mpFineMeshBoxCollection->GetLocalBoxes(boxForThisPoint);
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
            elem_index = mrFineMesh.GetContainingElementIndex(rPoint,
                                                              false,
                                                              test_element_indices,
                                                              true);
            weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

            mCounters[1]++;

        }
        catch(Exception& e)
        {
            if(safeMode)
            {
                // try the remaining elements
                try
                {
                    elem_index = mrFineMesh.GetContainingElementIndex(rPoint,
                                                                      false);
                    weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);
                    mCounters[2]++;

                }
                catch (Exception& e)
                {
                    // the point is not in ANY element, store the nearest element and corresponding weights
                    elem_index = mrFineMesh.GetNearestElementIndexFromTestElements(rPoint,test_element_indices);
                    weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

                    mNotInMesh.push_back(index);
                    mNotInMeshNearestElementWeights.push_back(weight);
                    mCounters[3]++;
                }
            }
            else
            {
                assert(test_element_indices.size()>0); // boxes probably too small if this fails
                
                // immediately assume it isn't in the rest of the mesh - this should be the 
                // case assuming the box width was chosen suitably.
                // Store the nearest element and corresponding weights
                elem_index = mrFineMesh.GetNearestElementIndexFromTestElements(rPoint,test_element_indices);
                weight = mrFineMesh.GetElement(elem_index)->CalculateInterpolationWeights(rPoint);

                mNotInMesh.push_back(index);
                mNotInMeshNearestElementWeights.push_back(weight);
                mCounters[3]++;
            }
        }
    }

    mFineMeshElementsAndWeights[index].ElementNum = elem_index;
    mFineMeshElementsAndWeights[index].Weights = weight;    
}


template<unsigned DIM>
void FineCoarseMeshPair<DIM>::PrintStatistics()
{
    assert(mNotInMesh.size()==mCounters[3]);
    assert(mNotInMesh.size()==mNotInMeshNearestElementWeights.size());
    std::cout << "\nFineCoarseMeshPair statistics:\n";
    std::cout << "\tNum points for which containing (fine) element was found, using box containing that point = " << mCounters[0] << "\n";
    std::cout << "\tNum points for which containing (fine) element was in local box = " << mCounters[1] << "\n";
    std::cout << "\tNum points for which containing (fine) element was in non-local element = " << mCounters[2] << "\n";
    std::cout << "\tNum points for which no containing element was found in fine mesh = " << mCounters[3] << "\n";
    if(mCounters[3]>0)
    {
        std::cout << "\tIndices and weights for points for which no containing element was found:\n";
        for(unsigned i=0; i<mNotInMesh.size(); i++)
        {
            std::cout << "\t\t" << mNotInMesh[i] << ", " << mNotInMeshNearestElementWeights[i] << "\n";
        }
    }
}



// todo #1409: use boxes as this is not coded in efficient way
template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeCoarseElementsForFineNodes()
{
    mCoarseElementsForFineNodes.resize(mrFineMesh.GetNumNodes());
    for(unsigned i=0; i<mCoarseElementsForFineNodes.size(); i++)
    {
        ChastePoint<DIM> point = mrFineMesh.GetNode(i)->GetPoint();

        try
        {
            mCoarseElementsForFineNodes[i] = mrCoarseMesh.GetContainingElementIndex(point);
        }
        catch(Exception& e)
        {
            mCoarseElementsForFineNodes[i] = mrCoarseMesh.GetNearestElementIndex(point);
        }
    }
}


// todo #1409: use boxes as this is not coded in efficient way
template<unsigned DIM>
void FineCoarseMeshPair<DIM>::ComputeCoarseElementsForFineElementCentroids()
{
    mCoarseElementsForFineElementCentroids.resize(mrFineMesh.GetNumElements());
    for(unsigned i=0; i<mrFineMesh.GetNumElements(); i++)
    {
        c_vector<double,DIM> point_cvec = mrFineMesh.GetElement(i)->CalculateCentroid();
        ChastePoint<DIM> point(point_cvec);
        try
        {
            mCoarseElementsForFineElementCentroids[i] = mrCoarseMesh.GetContainingElementIndex(point);
        }
        catch(Exception& e)
        {
            mCoarseElementsForFineElementCentroids[i] = mrCoarseMesh.GetNearestElementIndex(point);
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
