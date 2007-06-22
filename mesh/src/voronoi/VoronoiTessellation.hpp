#ifndef VORONOITESSELLATION_HPP_
#define VORONOITESSELLATION_HPP_

#include "UblasCustomFunctions.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "VoronoiCell.hpp"

#include <cmath>
#include <vector>

class VoronoiTessellation
{
  private:
    friend class TestVoronoiTessellation;
    
    ConformingTetrahedralMesh<3,3>& mrMesh;
    /**
     * Vertices correspond to elements of the mesh
     */
    std::vector< c_vector<double, 3>* > mVertices;
    /**
     * Faces corespond to edges of the mesh
     */
    std::vector< Face* > mFaces;
    /**
     * Cells correspond to nodes of the mesh
     */
    std::vector< VoronoiCell > mVoronoiCells;
  
  
    class VertexAndAngle
    {
    public:
        c_vector< double ,3 >* mpVertex;
        double mAngle;
        bool operator<(const VertexAndAngle& other) const
        {
            return mAngle < other.mAngle;
        }
    };
    
    void GenerateVerticesFromElementCircumcentres();
        
    double ReturnPolarAngle(double x, double y) const;
    


  public:
    
    
    /***
     * Constructor. Create a tesselation of the given mesh which must be Delauny
     * (see ConformingTetrahedralMesh::CheckVoronoi).
     */
    VoronoiTessellation(ConformingTetrahedralMesh<3,3>& rMesh);
    
    ~VoronoiTessellation();
    
    /***
     * Get a VoronoiCell.
     * 
     * @param index The index of the cell is the index of the corresponding node in the orginal mesh.
     * If the corresponding node was on the boundary, this will return a cell with no faces.
     */
    const VoronoiCell& rGetCell(unsigned index) const;

};

#endif /*VORONOITESSELLATION_HPP_*/
