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
#include "AbstractContractionModel.hpp"
#include "Kerchoffs2003ContractionModel.hpp"
#include "NonPhysiologicalContractionModel.hpp"


/**
 *  Explicit cardiac mechanics assembler for solving electromechanic problems where the
 *  contraction model is not stretch-rate-dependent (for those the implicit assembler is 
 *  needed). 
 *  
 *  The general explicit solution procedure is to do, each timestep:
 *  (0) [solve the electrics and interpolate Ca and voltage onto quad points
 *  (i) pass Ca and voltage to the contraction models
 *  (ii) pass the fibre stretch to the contraction models in case this is needed.
 *  (iii) integrate the contraction models in order to get the active tension
 *  (iv) solve for the deformation using this active tension. 
 */
template<unsigned DIM>
class ExplicitCardiacMechanicsAssembler : public AbstractCardiacMechanicsAssembler<DIM>
{
friend class TestExplicitCardiacMechanicsAssembler;

private:
    /**
     *  Vector of contraction model (pointers). One for each quadrature point.
     */
    std::vector<AbstractContractionModel*> mContractionModelSystems;
    
    /**
     *  Stored stretches (in fibre direction, at each quadrature point) from the 
     *  previous timestep, to pass to contraction models if needed.
     */
    std::vector<double> mStretches;

    /** This solver is an explicit solver (overloaded pure method) */
    bool IsImplicitSolver()
    {
        return false;
    }

    /**
     *  Get the active tension and other info at the given quadrature point. This is an explicit 
     *  assembler so just sets the active tension, it doesn't set the derivatives. It stores the
     *  stretch for the next timestep.
     * 
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex Quadrature point integrand currently being evaluated at in AssembleOnElement.
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension. 
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch. Unset in this explicit solver.
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate. Unset in this explicit solver.
     */   
    void GetActiveTensionAndTensionDerivs(double currentFibreStretch, 
                                          unsigned currentQuadPointGlobalIndex,
                                          bool assembleJacobian,
                                          double& rActiveTension,
                                          double& rDerivActiveTensionWrtLambda,
                                          double& rDerivActiveTensionWrtDLambdaDt)
    {
        // the active tensions have already been computed for each contraction model, so can 
        // return it straightaway..
        rActiveTension = mContractionModelSystems[currentQuadPointGlobalIndex]->GetActiveTension();

        // these are unset
        rDerivActiveTensionWrtLambda = 0.0;
        rDerivActiveTensionWrtDLambdaDt = 0.0;

        // store the value of given for this quad point, so that it can be used when computing 
        // the active tension at the next timestep
        mStretches[currentQuadPointGlobalIndex] = currentFibreStretch;
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
            case TEST1:
            {
                for(unsigned i=0; i<this->mTotalQuadPoints; i++)
                {
                    mContractionModelSystems.push_back(new NonPhysiologicalContractionModel(1));
                }
                break;
            }
            case KERCHOFFS2003: //stretch dependent, will this work with explicit??
            {
                for(unsigned i=0; i<this->mTotalQuadPoints; i++)
                {
                    Kerchoffs2003ContractionModel* p_model = new Kerchoffs2003ContractionModel();
                    mContractionModelSystems.push_back(p_model);
                }
                break;
            }
            default:
            {
                EXCEPTION("Unknown or stretch-rate-dependent contraction model");
            } 
        }
        
        mStretches.resize(this->mTotalQuadPoints);        
        assert(!(mContractionModelSystems[0]->IsStretchRateDependent()));
    }
    
    /**
     *  Destructor
     */
    virtual ~ExplicitCardiacMechanicsAssembler()
    {        
        for(unsigned i=0; i<mContractionModelSystems.size(); i++)
        {
            //// memory leak as this is commented out. But get glibc failure with it in... (EMTODO2)
            //delete mContractionModelSystems[i];
        }
    }        
        

    /**
     *  Set the intracellular Calcium concentrations and voltages at each quad point.
     * 
     *  This explicit solver (for contraction models which are NOT functions of stretch) can then
     *  integrate the contraction models to get the active tension, although this is done in Solve.
     * 
     *  @param rCalciumConcentrations Reference to a vector of intracellular calcium concentrations at each quadrature point
     *  @param rVoltages Reference to a vector of voltages at each quadrature point
     */

    void SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations, 
                              std::vector<double>& rVoltages)
    {
        assert(rCalciumConcentrations.size()==mContractionModelSystems.size());
        assert(rVoltages.size()==mContractionModelSystems.size());

        ContractionModelInputParameters input_parameters;
        
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

        // assemble the residual again so that mStretches is set (in GetActiveTensionAndTensionDerivs)
        // using the current deformation.
        this->AssembleSystem(true,false);
                
        // integrate contraction models
        for(unsigned i=0; i<mContractionModelSystems.size(); i++)
        {
            mContractionModelSystems[i]->SetStretchAndStretchRate(mStretches[i], 0.0);
            mContractionModelSystems[i]->RunDoNotUpdate(time, nextTime, odeTimestep);
            mContractionModelSystems[i]->UpdateStateVariables();
        }   
        
        // solve
        NonlinearElasticityAssembler<DIM>::Solve();
    }
};

#endif /*EXPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
