#ifndef ABSTRACTELASTICITYASSEMBLER_HPP_
#define ABSTRACTELASTICITYASSEMBLER_HPP_

#include "AbstractDealiiAssembler.hpp"
#include "DofVertexIterator.hpp"

template <unsigned DIM> 
class AbstractElasticityAssembler : public AbstractDealiiAssembler<DIM>
{
private:
    /*< Data structure containing the deformed position, by vertex index, in easily
     * accessable form. Only created if asked for */
    std::vector<Vector<double> > mDeformedPosition;

    /*< Data structure containing the undeformed position, by vertex index, in easily
     * accessable form. Only created if asked for */
    std::vector<Vector<double> > mUndeformedPosition;


public:
    /** Constructor
     *  
     *  Just takes in the mesh and passes it down to AbstractAssembler
     */
    AbstractElasticityAssembler(Triangulation<DIM>* pMesh)
        : AbstractDealiiAssembler<DIM>(pMesh)
    {
    }

    /**
     *  Get the deformed position. rGetDeformedPosition()[i][j] is the x_i value at node j
     */
    std::vector<Vector<double> >& rGetDeformedPosition()
    {
        mDeformedPosition.resize(DIM);
        for (unsigned i=0; i<DIM; i++)
        {
            mDeformedPosition[i].reinit(this->mpMesh->n_vertices());
        }
        
        DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<DIM> old_posn = vertex_iter.GetVertex();
            
            for (unsigned i=0; i<DIM; i++)
            {
                mDeformedPosition[i](vertex_index) =   old_posn(i)
                                                     + mCurrentSolution(vertex_iter.GetDof(i));
            }
            
            vertex_iter.Next();
        }
        
        return mDeformedPosition;
    }
    
    
    /**
     *  Get the undeformed position. rGetUndeformedPosition()[i][j] is the X_i value at node j
     *  Obviously this data is accessible from the mesh as well, this method is more useful
     *  in some situations. Note, this data structure is not set up unless this method 
     *  is called.
     * 
     *  Note we don't just calculate this once and store because the undeformed mesh will
     *  change if coarsening/refinement occurs
     */
    std::vector<Vector<double> >& rGetUndeformedPosition()
    {
        // initialise
        mUndeformedPosition.resize(DIM);
        for (unsigned i=0; i<DIM; i++)
        {
            mUndeformedPosition[i].reinit(this->mpMesh->n_vertices());
        }
        
        // populate
        TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<DIM> old_posn = vertex_iter.GetVertex();
            
            for (unsigned i=0; i<DIM; i++)
            {
                mUndeformedPosition[i](vertex_index) = vertex_iter.GetVertex()(i);
            }
            
            vertex_iter.Next();
        }
        
        return mUndeformedPosition;
    }
};
#endif /*ABSTRACTELASTICITYASSEMBLER_HPP_*/
