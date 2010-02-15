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
    FineCoarseMeshPair(TetrahedralMesh<DIM,DIM>& rFineMesh, QuadraticMesh<DIM>& rCoarseMesh);
    
    /**
     *  Destructor just deletes the box collection
     */
    ~FineCoarseMeshPair();
    
    /**
     *  Set up boxes on fine mesh. The elements contained in each box is stored, which makes
     *  finding the containing element for a given point much faster.
     *  This should be called before ComputeFineElementsAndWeightsForCoarseQuadPoints() etc
     * 
     *  @param boxWidth width to use for the boxes (which will be cubes). Note that a domain
     *    which is a touch larger than the smallest containing cuboid of the fine mesh is used.
     *    boxWidth defaults to a negative value, in which case a box width such that there are
     *    approximately 10 boxes in the x-direction.
     */
    void SetUpBoxesOnFineMesh(double boxWidth = -1);
    
    /**
     * Set up the containing (fine) elements and corresponding weights for all the  
     * quadrature points in the coarse mesh. Call GetElementsAndWeights() after calling this
     * with the index of the quad point (=the index of the quad point in a QuadraturePointsGroup=
     * the index if the quad points were listed by looping over all the element and then
     * looping over all the quad points). 
     * @rQuadRule The quadrature rule, used to determine the number of quadrature points per element.
     */
    void ComputeFineElementsAndWeightsForCoarseQuadPoints(GaussianQuadratureRule<DIM>& rQuadRule);
    
    /**
     *  Print the number of points for which the containing element was found quickly, the number
     *  for which the containing element was found slowly, and the number for which no containing 
     *  element was found (with the values of the weights for the latter).
     */
    void PrintStatistics();
    
    std::vector<ElementAndWeights<DIM> >& rGetElementsAndWeights()
    {
        return mElementsAndWeights;
    }
    
    /**
     *  Destroy the box collection - can be used to free memory once 
     *  ComputeFineElementsAndWeightsForCoarseQuadPoints has been called.
     */
    void DeleteBoxCollection();
};

#endif /*FINECOARSEMESHPAIR_HPP_*/
