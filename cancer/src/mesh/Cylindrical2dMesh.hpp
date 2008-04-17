#ifndef _CYLINDRICAL2DMESH_HPP_
#define _CYLINDRICAL2DMESH_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <math.h>

#include "ConformingTetrahedralMesh.hpp"
#include "NodeMap.hpp"
#include "TrianglesMeshWriter.hpp"

#include <boost/serialization/export.hpp>// at end of includes

class Cylindrical2dMesh : public ConformingTetrahedralMesh<2, 2>
{
    friend class TestCylindrical2dMesh;
private:
    
    /** The circumference of the cylinder */
    double mWidth;
    /** The top of the cylinder (y-coord) */
    double mTop;
    /** The bottom of the cylinder (y-coord) */
    double mBottom;
    
    /** Used during ReMesh **/
    // The left nodes which have been mirrored
    std::vector<unsigned> mLeftOriginals;
    // The image nodes relating to these left nodes (on right of mesh)
    std::vector<unsigned> mLeftImages;
    // The right nodes which have been mirrored
    std::vector<unsigned> mRightOriginals;
    // The image nodes relating to these right nodes (on left of mesh)
    std::vector<unsigned> mRightImages;
    // The indices of elements which straddle the left periodic boundary
    std::set<unsigned> mLeftPeriodicBoundaryElementIndices;
    // The indices of elements which straddle the right periodic boundary
    std::set<unsigned> mRightPeriodicBoundaryElementIndices;
    
    /** The indices of nodes on the top boundary */
    std::vector<unsigned > mTopHaloNodes;
    /** The indices of nodes on the bottom boundary */
    std::vector<unsigned > mBottomHaloNodes;
    
    void ReplaceImageWithRealNodeOnElement(Element<2,2>* pElement, std::vector<unsigned> &rImageNodes, std::vector<unsigned> &rOriginalNodes, unsigned nodeIndex ) ;
    
    // This bunch of functions should only ever be called by the public ReMesh
    void UpdateTopAndBottom();
    void CreateHaloNodes();
    void CreateMirrorNodes();
    void ReconstructCylindricalMesh();
    void DeleteHaloNodes();
    void CorrectNonPeriodicMesh();
    void GenerateVectorsOfElementsStraddlingPeriodicBoundaries();
    unsigned GetCorrespondingNodeIndex(unsigned nodeIndex);
    void UseTheseElementsToDecideMeshing(std::set<unsigned> mainSideElements);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<ConformingTetrahedralMesh<2,2> >(*this);
        archive & mWidth;
        archive & mTop;
        archive & mBottom;
    }
    
public:
    
    Cylindrical2dMesh(double width);
    Cylindrical2dMesh(double width, std::vector<Node<2> *> nodes);
         
    ~Cylindrical2dMesh()
    {
    }
        
    void ReMesh(NodeMap &map);
    c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2);
    void SetNode(unsigned index, ChastePoint<2> point, bool concreteMove);
    bool IsThisIndexInList(const unsigned& rNodeIndex, const std::vector<unsigned>& rListOfNodes);
    double GetWidth(const unsigned& rDimension) const;
    unsigned AddNode(Node<2> *pNewNode);
    
};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cylindrical2dMesh
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Cylindrical2dMesh * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const double width = t->GetWidth(0);
    ar << width;
}

/**
 * De-serialize constructor parameters and initialise Cylindrical2dMesh.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Cylindrical2dMesh * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    double width;
    ar >> width;
    // invoke inplace constructor to initialize instance
    ::new(t)Cylindrical2dMesh(width);
}
}
} // namespace ...

BOOST_CLASS_EXPORT(Cylindrical2dMesh)

#endif //_CYLINDRICAL2DMESH_HPP_
