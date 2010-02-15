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

#ifndef FINECOARSEMESHPAIR_HPP_
#define FINECOARSEMESHPAIR_HPP_

#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "BoxCollection.hpp"
#include "QuadraturePointsGroup.hpp"
#include "GaussianQuadratureRule.hpp"


/**
 *  At the beginning of a two mesh simulation we need to figure out and store
 *  which fine-mesh element each (coarse-mesh) quadrature point is in, and
 *  what the weight of that gauss point for that particular element is. This struct
 *  just contains this two pieces of data
 */
template<unsigned DIM>
struct ElementAndWeights
{
    unsigned ElementNum; /**< Which element*/
    c_vector<double,DIM+1> Weights; /**<Gauss weights for this element*/
};


/**
 *  Class for a pair of meshes, one fine, one coarse, which should cover the same domain (or very nearly match).
 *  This class is used to set up interpolation information from one mesh to the other
 */
template <unsigned DIM>
class FineCoarseMeshPair
{
friend class TestFineCoarseMeshPair;

private:
    /** Fine mesh */
    TetrahedralMesh<DIM,DIM>& mrFineMesh;
    /** Coarse mesh (usually be a quadratic mesh) */
    QuadraticMesh<DIM>& mrCoarseMesh;

    /** The min values of the nodes, for each dimension, in the fine mesh, for creating the boxes */    
    c_vector<double,DIM> mMinValuesFine;
    /** The max values of the nodes, for each dimension, in the fine mesh, for creating the boxes */    
    c_vector<double,DIM> mMaxValuesFine;

    /** Boxes on the fine mesh domain, for easier determination of containing element for a given point */
    BoxCollection<DIM>* mpFineMeshBoxCollection;
    /** The containing elements and corresponding weights in the fine mesh for the set of points given.
     *  The points may have been quadrature points in the coarse mesh, or nodes in coarse mesh, etc. 
     */
    std::vector<ElementAndWeights<DIM> > mElementsAndWeights;
    
    /** Indices of the points which were found to be outside the fine mesh */ 
    std::vector<unsigned> mNotInMesh;
    /** The corresponding weights, for the nearest elements, of the points which were found 
     *  to be outside the fine mesh */ 
    std::vector<c_vector<double,DIM+1> > mNotInMeshNearestElementWeights;
    
    /** 3 values, (0) number of points for which the containing element was found quickly (the element was
     *  in the same box as the point, (1) number of points for which the containing element was found 
     *  slowly (the element was not the same box as the point, (2) num points outside the fine mesh.
     *  Note mCounters[2] = mNotInMesh.size() = mNotInMeshNearestElementWeights.size();
     */
    std::vector<unsigned> mCounters;
    
////////////// only implement if going to be needed   
//////    /** In some simulations the coarse and fine meshes will turn out to be the same (the vertices of the
//////     *  coarse quadratic mesh will match the vertices of the fine mesh), we should figure out if this 
//////     *  is the case and do things differently if so
//////     */
//////    bool mIdenticalMeshes;

public:
    /** Constructor sets up domain size
     *  @param rFineMesh Fine mesh (reference)
     *  @param rCoarseMesh Coarse mesh (reference)
     */
    FineCoarseMeshPair(TetrahedralMesh<DIM,DIM>& rFineMesh, QuadraticMesh<DIM>& rCoarseMesh)
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
    
    /**
     *  Set up boxes on fine mesh. The elements contained in each box is stored, which makes
     *  finding the containing element for a given point much faster.
     *  This should be called before ComputeFineElementsAndWeightsForCoarseQuadPoints() etc
     *  @param boxWidth width to use for the boxes (which will be cubes). Note that a domain
     *  which is a touch larger than the smallest containing cuboid of the fine mesh is used.
     */
    void SetUpBoxesOnFineMesh(double boxWidth)
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
    
    /**
     * Set up the containing (fine) elements and corresponding weights for all the  
     * quadrature points in the coarse mesh. Call GetElementsAndWeights() after calling this
     * with the index of the quad point (=the index of the quad point in a QuadraturePointsGroup=
     * the index if the quad points were listed by looping over all the element and then
     * looping over all the quad points). 
     * @rQuadRule The quadrature rule, used to determine the number of quadrature points per element.
     */
    void ComputeFineElementsAndWeightsForCoarseQuadPoints(GaussianQuadratureRule<DIM>& rQuadRule)
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
    
    /**
     *  Print the number of points for which the containing element was found quickly, the number
     *  for which the containing element was found slowly, and the number for which no containing 
     *  element was found (with the values of the weights for the latter).
     */
    void PrintStatistics()
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
    
    std::vector<ElementAndWeights<DIM> >& rGetElementsAndWeights()
    {
        return mElementsAndWeights;
    }        
};

#endif /*FINECOARSEMESHPAIR_HPP_*/
