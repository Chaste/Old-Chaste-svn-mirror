#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include "AbstractCardiacCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"

/**
 * Class which returns cardiac cells (which may have been stimulated).
 * For use with MonodomainPde.
 * Needed because when running in parallel, the user should specify only
 * which cells are used in a mesh, not on which process they exist.
 * A concrete implementation will be required, possibly for each simulation.
 */
 
template<int SPACE_DIM> 
class AbstractCardiacCellFactory
{
protected:
    double mTimeStep;
    InitialStimulus* mpZeroStimulus;
    AbstractIvpOdeSolver* mpSolver;
 
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

public:
    virtual AbstractCardiacCell* CreateCardiacCellForNode(int)=0;
    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed, int lo, int hi) {}

    virtual int GetNumberOfNodes()
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
        mpZeroStimulus = new InitialStimulus(0,0,0);
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

