/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef EXPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define EXPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include "AbstractCardiacMechanicsAssembler.hpp"
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractOdeBasedContractionModel.hpp"

#include "Nash2004ContractionModel.hpp"
#include "Kerchoffs2003ContractionModel.hpp"


/**
 */
template<unsigned DIM>
class ExplicitCardiacMechanicsAssembler : public AbstractCardiacMechanicsAssembler<DIM>
{
friend class TestExplicitCardiacMechanicsAssembler;

private:
    std::vector<AbstractOdeBasedContractionModel*> mContractionModelSystems;
    std::vector<double> mActiveTensions;

    /** This solver is an explicit solver (overloaded pure method) */
    bool IsImplicitSolver()
    {
        return false;
    }

    void GetActiveTensionAndTensionDerivs(c_matrix<double,DIM,DIM>& C, 
                                          unsigned currentQuadPointGlobalIndex,
                                          bool assembleJacobian,
                                          double& rActiveTension,
                                          double& rDerivActiveTensionWrtLambda,
                                          double& rDerivActiveTensionWrtDLambdaDt,
                                          double& rLambda)
    {
        rActiveTension = mActiveTensions[currentQuadPointGlobalIndex];
    }

public:
    /**
     * Constructor
     *
     * @param pQuadMesh A pointer to the mesh.
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     * @param rFixedNodes The fixed nodes
     * @param pMaterialLaw The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    ExplicitCardiacMechanicsAssembler(ContractionModel contractionModel,
                                      QuadraticMesh<DIM>* pQuadMesh,
                                      std::string outputDirectory,
                                      std::vector<unsigned>& rFixedNodes,
                                      AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw = NULL)
        : AbstractCardiacMechanicsAssembler<DIM>(pQuadMesh,
                                                 outputDirectory,
                                                 rFixedNodes,
                                                 pMaterialLaw)
    {
        switch(contractionModel)
        {
            case NASH2004:
            {
                mContractionModelSystems.resize(this->mTotalQuadPoints, new Nash2004ContractionModel());
                break;
            }
            case KERCHOFFS2003: //stretch dependent, will this work with explicit??
            {
                mContractionModelSystems.resize(this->mTotalQuadPoints, new Nash2004ContractionModel());
                break;
            }
            default:
            {
                EXCEPTION("Unknown or stretch-rate-dependent contraction model");
            } 
        }
        
        assert(!(mContractionModelSystems[0]->IsStretchRateDependent()));
        
        // initialise stores
        mActiveTensions.resize(this->mTotalQuadPoints, 0.0);
    }
    /**
     *  Destructor
     */
    virtual ~ExplicitCardiacMechanicsAssembler()
    {
        for(unsigned i=0; i<mContractionModelSystems.size(); i++)
        {
            delete mContractionModelSystems[i];
        }
    }        
        

    /**
     *  Set the intracellular Calcium concentrations (note: in an explicit algorithm we
     *  would set the active tension as the forcing quantity; the implicit algorithm
     *  takes in the Calcium concentration and solves for the active tension implicitly
     *  together with the mechanics).
     * 
     *  @param caI the intracellular calcium concentrations
     */
    void SetIntracellularCalciumConcentrations(std::vector<double>& caI)
    {
        assert(caI.size()==mContractionModelSystems.size());
        assert(caI.size()==mActiveTensions.size());

        ContractionModelInputParameters input_parameters;
        input_parameters.Voltage = DOUBLE_UNSET;
        input_parameters.Time = DOUBLE_UNSET;
        
        for(unsigned i=0; i<mContractionModelSystems.size(); i++)
        {
            input_parameters.IntracellularCalciumConcentration = caI[i];
            mContractionModelSystems[i]->SetInputParameters(input_parameters);
        }
    }
    
    /**
     *  Solve for the deformation using quasi-static nonlinear elasticity.
     *  (not dynamic nonlinear elasticity, despite the times taken in - just ONE
     *  deformation is solved for. The cell models are integrated explicitly
     *  over the time range using the ODE timestep provided then the active tension
     *  used to solve for the deformation
     * 
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    void Solve(double time, double nextTime, double odeTimestep)
    {
        assert(time < nextTime);
        this->mCurrentTime = time;
        this->mNextTime = nextTime;
        this->mOdeTimestep = odeTimestep;        
                
        EulerIvpOdeSolver euler_solver;
        for(unsigned i=0; i<mContractionModelSystems.size(); i++)
        {
            euler_solver.SolveAndUpdateStateVariable(mContractionModelSystems[i], time, nextTime, odeTimestep);
            mActiveTensions[i] = mContractionModelSystems[i]->GetActiveTension();
        }   
        
        // solve
        NonlinearElasticityAssembler<DIM>::Solve();
    }
};

#endif /*EXPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
