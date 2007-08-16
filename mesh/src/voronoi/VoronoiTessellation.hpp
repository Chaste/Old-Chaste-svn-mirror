#ifndef VORONOITESSELLATION_HPP_
#define VORONOITESSELLATION_HPP_

#include "UblasCustomFunctions.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "VoronoiCell.hpp"

#include <cmath>
#include <vector>

template<unsigned DIM>
class VoronoiTessellation
{
  private:
    friend class TestVoronoiTessellation;
    friend class InventorVoronoiWriter;
    
    ConformingTetrahedralMesh<DIM,DIM>& mrMesh;
    /**
     * Vertices correspond to elements of the mesh
     */
    std::vector< c_vector<double,DIM>* > mVertices;
    /**
     * Faces corespond to edges of the mesh
     */
    std::vector< Face<DIM>* > mFaces;
    /**
     * Cells correspond to nodes of the mesh
     */
    std::vector< VoronoiCell > mVoronoiCells;
  
  
    class VertexAndAngle
    {
    public:
        c_vector< double ,DIM >* mpVertex;
        double mAngle;
        bool operator<(const VertexAndAngle& other) const
        {
            return mAngle < other.mAngle;
        }
    };
    
    void GenerateVerticesFromElementCircumcentres();
    
    /**
     * @param x x-coordinate
     * @param y y-coordinate
     * @return Polar angle in interval (-PI,PI]
     */    
    double ReturnPolarAngle(double x, double y) const;
    
    void Initialise(ConformingTetrahedralMesh<2,2>& rMesh);
    void Initialise(ConformingTetrahedralMesh<3,3>& rMesh);
    

  public:
    
    
    /***
     * Constructor. Create a tesselation of the given mesh which must be Delauny
     * (see ConformingTetrahedralMesh::CheckVoronoi).
     */
    VoronoiTessellation(ConformingTetrahedralMesh<DIM,DIM>& rMesh);
    
    ~VoronoiTessellation();
    
    /***
     * Get a VoronoiCell.
     * 
     * @param index The index of the cell is the index of the corresponding node in the orginal mesh.
     * If the corresponding node was on the boundary, this will return a cell with no faces.
     */
    const VoronoiCell& rGetCell(unsigned index) const;
    const Face<DIM>* GetFace(unsigned index) const;
    
};

#endif /*VORONOITESSELLATION_HPP_*/
