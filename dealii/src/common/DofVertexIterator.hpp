#ifndef DOFVERTEXITERATOR_HPP_
#define DOFVERTEXITERATOR_HPP_

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <cassert>

/** 
 *  DofVertexIterator
 * 
 *  See also TriangulationVertexIterator
 * 
 *  Looping over nodes (vertices) in deal.II meshes is a pain, you have to loop over
 *  elements (cells), then loop over each vertex of that cell, checking you haven't
 *  already done that vertex. This class takes care of this for you. 
 *  
 *  The difference between this class and TriangulationVertexIterator is that the cell
 *  is obtained through a DoFHandler<DIM>::active_cell_iterator rather than a 
 *  Triangulation<DIM>::active_cell_iterator, which means the dofs for the vertex
 *  can be obtained.
 * 
 *  Note: this class doesn't act like a typical std::iterator class, in that the
 *  result doesn't point to a vertex, you don't call ++ etc. (although in the 
 *  future ++ might be implemented to call Next()
 * 
 *  Usage:
 * 
 *  DofVertexIterator<2> iter(&mesh,&dof);
 *  while(!iter.ReachedEnd())
 *  {
 *      Point<2> vertex = iter.GetVertex(); //etc
 * 
 *      double value_at_node = solutions( iter.GetDof(0) ); // here solutions is a Vector<double>
 * 
 *      iter.Next();
 *  }
 * 
 *  The vertex number, current cell, local number of the vertex in that cell can
 *  also be obtained. 
 * 
 *  TODO: extract commonality between this and TriangulationVertexIterator. (although 
 *  the code  is almost identical mpCell are different types in the two classes)
 */
template<unsigned DIM>
class DofVertexIterator
{
private :
    Triangulation<DIM>* mpMesh;
    DoFHandler<DIM>* mpDofHandler;
    std::vector<bool> mVertexTouched;
    typename DoFHandler<DIM>::active_cell_iterator mpCurrentCell;

    unsigned mCurrentVertexIndex;
    bool mReachedEnd;
    
    void NextNode()
    {
        if(mCurrentVertexIndex < GeometryInfo<DIM>::vertices_per_cell-1)
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
    DofVertexIterator(Triangulation<DIM>* pMesh, DoFHandler<DIM>* pDofHandler)
       : mpMesh(pMesh),
         mpDofHandler(pDofHandler),
         mVertexTouched(mpMesh->n_vertices(),false),
         mpCurrentCell(mpDofHandler->begin_active())
    {
        assert(pMesh && pDofHandler);

        mCurrentVertexIndex = 0;
        
        mReachedEnd = (mpCurrentCell==mpDofHandler->end());
    }
    
    /**
     *  Move to the next vertex
     */
    void Next()
    {
        bool found=false;
        while( (found==false) && (!mReachedEnd) )
        {
            NextNode();

            mReachedEnd = (mpCurrentCell==mpDofHandler->end());
            if( !mReachedEnd && !mVertexTouched[mpCurrentCell->vertex_index(mCurrentVertexIndex)] )
            {
                found = true;
            }               
        }
        
        if(!mReachedEnd)
        {
            mVertexTouched[mpCurrentCell->vertex_index(mCurrentVertexIndex)] = true;
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
    typename DoFHandler<DIM>::active_cell_iterator GetCell()
    {
        assert(!mReachedEnd);
        return mpCurrentCell;
    }        
     
    /** 
     * Get the Dof associated with the ith unknown at the current vertex
     */
    unsigned GetDof(unsigned i)
    {      
        return mpCurrentCell->vertex_dof_index(mCurrentVertexIndex,i);
    }
       
    /** 
     *  Reset to the first vertex so the iterator can be used again
     */    
    void Reset()
    {
        for(unsigned i=0; i<mpMesh->n_vertices(); i++)
        {
            mVertexTouched[i] = false;
        }
        mpCurrentCell = mpDofHandler->begin_active();
        mCurrentVertexIndex = 0;
        mReachedEnd = (mpCurrentCell==mpDofHandler->end());
    }
};

#endif /*DOFVERTEXITERATOR_HPP_*/
