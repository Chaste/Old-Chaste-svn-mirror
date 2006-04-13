#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include "AbstractCardiacCell.hpp"
#include "ConformingTetrahedralMesh.hpp"

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
    virtual int GetNumberOfNodes()=0;
    virtual ~AbstractCardiacCellFactory() {}
    void SetMesh( ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
    {
        mpMesh = pMesh;
    }
    
};

#endif /*ABSTRACTCARDIACCELLFACTORY_HPP_*/

