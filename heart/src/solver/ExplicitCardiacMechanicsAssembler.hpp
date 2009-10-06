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
 *  Explicit cardiac mechanics assembler for solving electromechanic problems where the
 *  contraction model is not stretch-rate-dependent (for those the implicit assembler is 
 *  needed). 
 *  
 *  The general explicit solution procedure is to do, each timestep:
 *  (0) [solve the electrics and interpolate Ca and voltage onto quad points
 *  (i) pass Ca and voltage to the contraction models
 *  (ii) integrate the contraction models in order to get the active tension
 *  (iii) solve for the deformation using this active tension. 
 */
template<unsigned DIM>
class ExplicitCardiacMechanicsAssembler : public AbstractCardiacMechanicsAssembler<DIM>
{
friend class TestExplicitCardiacMechanicsAssembler;

private:
    /**
     *  Vector of contraction model (pointers). One for each quadrature point.
     */
    std::vector<AbstractOdeBasedContractionModel*> mContractionModelSystems;

    /** This solver is an explicit solver (overloaded pure method) */
    bool IsImplicitSolver()
    {
        return false;
    }

    /**
     *  Get the active tension and other info at the given quadrature point. This is an explicit 
     *  assembler so just sets the active tension, it doesn't set the derivatives or the stretch.
     * 
     *  @param C Green-deformation tension. Unused.
     *  @param currentQuadPointGlobalIndex quadrature point integrand currently being evaluated at in AssembleOnElement.
     *  @param rActiveTension The returned active tension. 
     *  @param rActiveTension The returned dT_dLam, derivative of active tension wrt stretch. Unset.
     *  @param rActiveTension The returned dT_dLamDot, derivative of active tension wrt stretch rate. Unset.
     *  @param rLambda The stretch (computed from C) as AssembleOnElement needs to use this too. Unset.
     */   
    void GetActiveTensionAndTensionDerivs(c_matrix<double,DIM,DIM>& C, 
                                          unsigned currentQuadPointGlobalIndex,
                                          bool assembleJacobian,
                                          double& rActiveTension,
                                          double& rDerivActiveTensionWrtLambda,
                                          double& rDerivActiveTensionWrtDLambdaDt,
                                          double& rLambda)
    {
        rActiveTension = mContractionModelSystems[currentQuadPointGlobalIndex]->GetActiveTension();
    }

public:
    /**
     * Constructor
     *
     * @param contractionModel The contraction model.
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
                mContractionModelSystems.resize(this->mTotalQuadPoints, new Kerchoffs2003ContractionModel());
                break;
            }
            default:
            {
                EXCEPTION("Unknown or stretch-rate-dependent contraction model");
            } 
        }
        
        assert(!(mContractionModelSystems[0]->IsStretchRateDependent()));
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
     *  Set the intracellular Calcium concentrations and voltages at each quad point, and the current time.
     * 
     *  This explicit solver (for contraction models which are NOT functions of stretch) can then
     *  integrate the contraction models to get the active tension, although this is done in Solve.
     * 
     *  @param rCalciumConcentrations Reference to a vector of intracellular calcium concentrations at each quadrature point
     *  @param rVoltages Reference to a vector of voltages at each quadrature point
     *  @param time Current time
     */

    void SetCalciumVoltageAndTime(std::vector<double>& rCalciumConcentrations, 
                                  std::vector<double>& rVoltages,
                                  double time)
    {
        assert(rCalciumConcentrations.size()==mContractionModelSystems.size());
        assert(rVoltages.size()==mContractionModelSystems.size());

        ContractionModelInputParameters input_parameters;

        input_parameters.time = time;
        
        for(unsigned i=0; i<mContractionModelSystems.size(); i++)
        {
            input_parameters.intracellularCalciumConcentration = rCalciumConcentrations[i];
            input_parameters.voltage = rVoltages[i];
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
        }   
        
        // solve
        NonlinearElasticityAssembler<DIM>::Solve();
    }
};

#endif /*EXPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
