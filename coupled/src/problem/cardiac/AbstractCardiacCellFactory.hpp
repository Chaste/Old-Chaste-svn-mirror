#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include "AbstractCardiacCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "ZeroStimulus.hpp"
#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"

/**
 * Class which returns cardiac cells.
 * For use with MonodomainPde and BidomainPdes.
 *
 * The user should implement their own concrete class, in particular implementing
 * CreateCardiacCellForNode(unsigned), which should return the cell corresponding to a
 * given node. The user should also implement GetNumberOfCells() if this isn't equal
 * to the number of nodes. FinaliseCellCreation() can be used to (eg) add stimuli to
 * certain cells after they have been created.
 *
 * This class saves the user having to create cells in parallel, that work is done
 * by the pde instead.
 */

template<unsigned SPACE_DIM>
class AbstractCardiacCellFactory
{
protected:
    double mTimeStep;
    ZeroStimulus* mpZeroStimulus;
    AbstractIvpOdeSolver* mpSolver;
    
    /** the mesh is automatically set in MonodomainProblem and BidomainProblem */
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;
    
public:
    virtual AbstractCardiacCell* CreateCardiacCellForNode(unsigned)=0;
    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {}
    
    virtual unsigned GetNumberOfCells()
    {
        assert(mpMesh != NULL);
        return mpMesh->GetNumNodes();
    }
    
    AbstractCardiacCellFactory(double timeStep,
                               AbstractIvpOdeSolver* pSolver = new EulerIvpOdeSolver)
    {
        mTimeStep = timeStep;
        mpMesh = NULL;
        mpSolver = pSolver;
        mpZeroStimulus = new ZeroStimulus;
    }
    virtual ~AbstractCardiacCellFactory()
    {
        delete mpSolver;
        delete mpZeroStimulus;
    }
    void SetMesh(ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
    {
        mpMesh = pMesh;
    }
    
};

#endif /*ABSTRACTCARDIACCELLFACTORY_HPP_*/

