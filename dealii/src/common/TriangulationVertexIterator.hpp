#ifndef TRIANGULATIONVERTEXITERATOR_HPP_
#define TRIANGULATIONVERTEXITERATOR_HPP_

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <cassert>


/**
 *  TriangulationVertexIterator
 *
 *  Looping over nodes (vertices) in deal.II meshes is a pain, you have to loop over
 *  elements (cells), then loop over each vertex of that cell, checking you haven't
 *  already done that vertex. This class takes care of this for you.
 *
 *  Note: this class doesn't act like a typical std::iterator class, in that the
 *  result doesn't point to a vertex, you don't call ++ etc. (although in the
 *  future ++ might be implemented to call Next()
 *
 *  See also DofVertexIterator
 *
 *  Usage:
 *
 *  TriangulationVertexIterator<2> iter(&mesh)
 *  while(!iter.ReachedEnd())
 *  {
 *      Point<2> vertex = iter.GetVertex(); //etc
 *
 *      iter.Next();
 *  }
 *
 *  The vertex number, current cell, local number of the vertex in that cell can
 *  also be obtained.
 * 
 *  Note: vertices are reached by looping over elements (internally, in this class), 
 *  therefore this class only reaches ACTIVE vertices (ie ones that are part of an
 *  element). (Unactive vertices arise if the mesh is coarsened)
 *
 *  TODO: perhaps extract commonality between this and DofVertexIterator. (although 
 *  the code is almost identical mpCell are different types in the two classes)
 */
template<unsigned DIM>
class TriangulationVertexIterator
{
private :
    Triangulation<DIM>* mpMesh;
    std::vector<bool> mVertexTouched;
    typename Triangulation<DIM>::active_cell_iterator mpCurrentCell;
    
    unsigned mCurrentVertexIndex;
    bool mReachedEnd;
    
    void NextNode()
    {
        if (mCurrentVertexIndex < GeometryInfo<DIM>::vertices_per_cell-1)
        {
            mCurrentVertexIndex++;
        }
        else
        {
            mCurrentVertexIndex = 0;
            mpCurrentCell++;
        }
    }
    
public :
    TriangulationVertexIterator(Triangulation<DIM>* pMesh)
            : mpMesh(pMesh),
            mVertexTouched(mpMesh->n_vertices(),false),
            mpCurrentCell(mpMesh->begin_active())
    {
        assert(pMesh);
        
        mCurrentVertexIndex = 0;
        
        mReachedEnd = (mpCurrentCell==mpMesh->end());

        // set the current node as having been touched
        if(!mReachedEnd)
        {
            mVertexTouched[GetVertexGlobalIndex()] = true;
        }
   }
    
    /**
      *  Move to the next vertex
      */
    void Next()
    {
        bool found=false;
        while ( (found==false) && (!mReachedEnd) )
        {
            NextNode();
            
            mReachedEnd = (mpCurrentCell==mpMesh->end());
            
            if ( !mReachedEnd && !mVertexTouched[mpCurrentCell->vertex_index(mCurrentVertexIndex)] )
            {
                found = true;
                mVertexTouched[GetVertexGlobalIndex()] = true;
            }
        }
    }
    
    
    /**
     *  The method returns true if the last vertex has been passed
     */
    bool ReachedEnd()
    {
        return mReachedEnd;
    }
    
    /**
     *  Get the position of the vertex
     */
    Point<DIM>& GetVertex()
    {
        assert(!mReachedEnd);
        return mpCurrentCell->vertex(mCurrentVertexIndex);
    }
    
    /**
     *  Get the global vertex index
     */
    unsigned GetVertexGlobalIndex()
    {
        assert(!mReachedEnd);
        return mpCurrentCell->vertex_index(mCurrentVertexIndex);
    }
    
    /**
     *  Get the index of the current vertex in the current cell
     *  To be used with GetCell()
     */
    unsigned GetLocalVertexIndexForCell()
    {
        assert(!mReachedEnd);
        return mCurrentVertexIndex;
    }
    
    /**
     *  Get the current cell. Together with GetLocalVertexIndexForCell() this
     *  specifies the vertex. Needed for calling other methods on the cell. 
     */
    typename Triangulation<DIM>::active_cell_iterator GetCell()
    {
        assert(!mReachedEnd);
        return mpCurrentCell;
    }
    
    /**
     *  Reset to the first vertex so the iterator can be used again
     */
    void Reset()
    {
        // resize mVertexTouched in case the mesh has been refined..
        if(mVertexTouched.size()!=mpMesh->n_vertices())
        {
            mVertexTouched.resize(mpMesh->n_vertices());
        }

        for (unsigned i=0; i<mpMesh->n_vertices(); i++)
        {
            mVertexTouched[i] = false;
        }
        mpCurrentCell = mpMesh->begin_active();
        mCurrentVertexIndex = 0;
        mReachedEnd = (mpCurrentCell==mpMesh->end());
        
        // set the current node as having been touched
        if(!mReachedEnd)
        {
            mVertexTouched[GetVertexGlobalIndex()] = true;
        }
    }
};


#endif /*TRIANGULATIONVERTEXITERATOR_HPP_*/
