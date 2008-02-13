#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

#include <boost/numeric/ublas/matrix.hpp>

#include "MonodomainDg0Assembler.hpp"
#include "MonodomainPde.hpp"
#include "AbstractCardiacProblem.hpp"


/**
 * Class which specifies and solves a monodomain problem.
 */
template<unsigned SPACE_DIM>
class MonodomainProblem : public AbstractCardiacProblem<SPACE_DIM, 1>
{
private:
    MonodomainPde<SPACE_DIM>* mpMonodomainPde;    

public:
    AbstractCardiacPde<SPACE_DIM>* CreateCardiacPde()
    {
        mpMonodomainPde = new MonodomainPde<SPACE_DIM>(this->mpCellFactory);

        this->mIntracellullarConductivityTensors.Init();                
        mpMonodomainPde->SetIntracellularConductivityTensors( & this->mIntracellullarConductivityTensors );

        return mpMonodomainPde;
    }
    
    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, 1>* CreateAssembler()
    {
        MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM>* p_assembler 
          = new MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM>(this->mpMesh, mpMonodomainPde, 2);
        p_assembler->SetBoundaryConditionsContainer(this->mpBoundaryConditionsContainer);
        if (this->mUseLinearSolverAbsoluteTolerance)
        {
            p_assembler->SetLinearSolverAbsoluteTolerance(this->mLinearSolverTolerance);       
        }
        else
        {
            p_assembler->SetLinearSolverRelativeTolerance(this->mLinearSolverTolerance);    
        }
        
        return p_assembler;
    }

public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : AbstractCardiacProblem<SPACE_DIM, 1>(pCellFactory),
            mpMonodomainPde(NULL)
    {
    }
    
    /**
     * Destructor
     */
    ~MonodomainProblem()
    {
    }
    
    MonodomainPde<SPACE_DIM> * GetMonodomainPde()
    {
        assert(mpMonodomainPde != NULL);
        return mpMonodomainPde;
    }
    
    /**
     *  Print out time and max/min voltage values at current time.
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(this->mVoltage);
        double v_max = -1e5, v_min = 1e5;
        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            if ( voltage_replicated[i] > v_max)
            {
                v_max = voltage_replicated[i];
            }
            if ( voltage_replicated[i] < v_min)
            {
                v_min = voltage_replicated[i];
            }
        }
        std::cout << " max/min V = " <<  v_max << " " <<   v_min << "\n" << std::flush;
    }
};

#endif /*MONODOMAINPROBLEM_HPP_*/
