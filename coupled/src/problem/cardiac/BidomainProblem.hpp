#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_



#include "BidomainDg0Assembler.hpp"
#include "BidomainPde.hpp"
#include "AbstractCardiacProblem.hpp"


/**
 * Class which specifies and solves a bidomain problem.
 *
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the
 * extracellular potential at node j.
 */
template<unsigned SPACE_DIM>
class BidomainProblem : public AbstractCardiacProblem<SPACE_DIM, 2>
{
private:    
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    std::vector<unsigned> mFixedExtracellularPotentialNodes; /** nodes at which the extracellular voltage is fixed to zero (replicated) */    
    unsigned mExtracelluarColumnId;
    
protected:
    AbstractCardiacPde<SPACE_DIM> *CreateCardiacPde()
    {
        mpBidomainPde = new BidomainPde<SPACE_DIM>(this->mpCellFactory);
        return mpBidomainPde;
    }
    
    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, 2>* CreateAssembler()
    {
        BidomainDg0Assembler<SPACE_DIM,SPACE_DIM>* p_bidomain_assembler
            = new BidomainDg0Assembler<SPACE_DIM,SPACE_DIM>(this->mpMesh, mpBidomainPde,
                                                            2, this->mLinearSolverRelativeTolerance);
        try
        {
            p_bidomain_assembler->SetFixedExtracellularPotentialNodes(mFixedExtracellularPotentialNodes);
        }
        catch (const Exception& e)
        {
            delete p_bidomain_assembler;
            throw e;
        }
        
        return p_bidomain_assembler;
    }
    
public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : AbstractCardiacProblem<SPACE_DIM, 2>(pCellFactory),
            mpBidomainPde(NULL)
    {
        mFixedExtracellularPotentialNodes.resize(0);
    }
    
    /**
     * Destructor
     */
    ~BidomainProblem()
    {
    }
    
    
    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to 
     *  zero. This does not necessarily have to be called. If it is not, phi_e 
     *  is only defined up to a constant.
     * 
     *  @param the nodes to be fixed.
     * 
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes)
    {        
        mFixedExtracellularPotentialNodes.resize(nodes.size());
        for (unsigned i=0; i<nodes.size(); i++)
        {
            // the assembler checks that the nodes[i] is less than
            // the number of nodes in the mesh so this is not done here
            mFixedExtracellularPotentialNodes[i] = nodes[i];
        }
    }
    
    /**
     *  Get the pde. Can only be called after Initialise()
     */
    BidomainPde<SPACE_DIM>* GetBidomainPde()
    {
        assert(mpBidomainPde!=NULL);
        return mpBidomainPde;
    }
    
    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(this->mVoltage);
        double v_max = -1e5, v_min = 1e5, phi_max = -1e5, phi_min = 1e5;

        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            if ( voltage_replicated[2*i] > v_max)
            {
                v_max = voltage_replicated[2*i];
            }
            if ( voltage_replicated[2*i] < v_min)
            {
                v_min = voltage_replicated[2*i];
            }
            if ( voltage_replicated[2*i+1] > phi_max)
            {
                phi_max = voltage_replicated[2*i+1];
            }
            if ( voltage_replicated[2*i+1] < phi_min)
            {
                phi_min = voltage_replicated[2*i+1];
            }
        }
        std::cout << " max/min V, phi_e = " << v_max << " " << v_min << " " << phi_max << " " << phi_min << "\n" << std::flush;
    }
    
    virtual void DefineWriterColumns()
    {
        AbstractCardiacProblem<SPACE_DIM,2>::DefineWriterColumns();
        mExtracelluarColumnId = this->mpWriter->DefineVariable("Phi_e","mV");
    }
    
    virtual void WriteOneStep(double time, Vec voltageVec)
    {
        AbstractCardiacProblem<SPACE_DIM,2>::WriteOneStep(time, voltageVec);
        DistributedVector::Stripe extracellular(voltageVec, 1);
        this->mpWriter->PutVectorStripe(mExtracelluarColumnId, extracellular);
    }
    
};


#endif /*BIDOMAINPROBLEM_HPP_*/
