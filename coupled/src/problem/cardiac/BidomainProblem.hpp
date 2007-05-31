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
    double mLinearSolverRelativeTolerance;
    
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
        mLinearSolverRelativeTolerance=1e-6;
    }
    
    /**
     * Destructor
     */
    ~BidomainProblem()
    {
    }
    
    AbstractCardiacPde<SPACE_DIM>* CreatePde()
    {
        mpBidomainPde = new BidomainPde<SPACE_DIM>( this->mpCellFactory );
        return mpBidomainPde;
    }
    
    AbstractLinearDynamicProblemAssembler<SPACE_DIM, SPACE_DIM, 2>* CreateAssembler()
    {
         BidomainDg0Assembler<SPACE_DIM,SPACE_DIM>* p_bidomain_assembler = 
                 new BidomainDg0Assembler<SPACE_DIM,SPACE_DIM>(&this->mMesh, mpBidomainPde);
         p_bidomain_assembler->SetFixedExtracellularPotentialNodes(mFixedExtracellularPotentialNodes);
         return p_bidomain_assembler;
    }
    
    void SetLinearSolverRelativeTolerance(const double &rRelTol)
    {
        mLinearSolverRelativeTolerance = rRelTol;
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

        for (unsigned i=0; i<this->mMesh.GetNumNodes(); i++)
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
        std::cout << " max/min V, phi_e = "
        << v_max << " "
        << v_min << " "
        << phi_max << " "
        << phi_min << "\n" << std::flush;
    }
    
    std::string ColumnName()
    {
        return "Vm_And_Phi_e";
    }
};


#endif /*BIDOMAINPROBLEM_HPP_*/
