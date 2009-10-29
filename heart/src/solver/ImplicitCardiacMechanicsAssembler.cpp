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

#include "ImplicitCardiacMechanicsAssembler.hpp"

template<unsigned DIM>
ImplicitCardiacMechanicsAssembler<DIM>::ImplicitCardiacMechanicsAssembler(
                                  QuadraticMesh<DIM>* pQuadMesh,
                                  std::string outputDirectory,
                                  std::vector<unsigned>& rFixedNodes,
                                  AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
    : AbstractCardiacMechanicsAssembler<DIM>(pQuadMesh,
                                             outputDirectory,
                                             rFixedNodes,
                                             pMaterialLaw)
{
    // initialise stores
    mLambda.resize(this->mTotalQuadPoints, 1.0);
    mLambdaLastTimeStep.resize(this->mTotalQuadPoints, 1.0);
    mCellMechSystems.resize(this->mTotalQuadPoints);
}

template<unsigned DIM>
ImplicitCardiacMechanicsAssembler<DIM>::~ImplicitCardiacMechanicsAssembler()
{
}


template<unsigned DIM>
void ImplicitCardiacMechanicsAssembler<DIM>::SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations, 
                                                                  std::vector<double>& rVoltages)
                                        
{
    assert(rCalciumConcentrations.size() == this->mTotalQuadPoints);
    assert(rVoltages.size() == this->mTotalQuadPoints);

    ContractionModelInputParameters input_parameters;
    
    for(unsigned i=0; i<rCalciumConcentrations.size(); i++)
    {
        input_parameters.intracellularCalciumConcentration = rCalciumConcentrations[i];
        input_parameters.voltage = rVoltages[i];
        
        mCellMechSystems[i].SetInputParameters(input_parameters);
    }
}

template<unsigned DIM>
std::vector<double>& ImplicitCardiacMechanicsAssembler<DIM>::rGetLambda()
{
    return mLambda;
}


template<unsigned DIM>
void ImplicitCardiacMechanicsAssembler<DIM>::Solve(double time, double nextTime, double odeTimestep)
{
    // set the times, which are used in AssembleOnElement
    assert(time < nextTime);
    this->mCurrentTime = time;
    this->mNextTime = nextTime;
    this->mOdeTimestep = odeTimestep;

    // solve
    NonlinearElasticityAssembler<DIM>::Solve();

    // assemble residual again (to solve the cell models implicitly again
    // using the correct value of the deformation x (in case this wasn't the
    // last thing that was done
    this->AssembleSystem(true,false);

    // now update state variables, and set lambda at last timestep. Note
    // lambda was set in AssembleOnElement
    for(unsigned i=0; i<mCellMechSystems.size(); i++)
    {
         mCellMechSystems[i].UpdateStateVariables();
         mLambdaLastTimeStep[i] = mLambda[i];
    }
}



template<unsigned DIM>
void ImplicitCardiacMechanicsAssembler<DIM>::GetActiveTensionAndTensionDerivs(double currentFibreStretch, 
                                                                              unsigned currentQuadPointGlobalIndex,
                                                                              bool assembleJacobian,
                                                                              double& rActiveTension,
                                                                              double& rDerivActiveTensionWrtLambda,
                                                                              double& rDerivActiveTensionWrtDLambdaDt)
{
    // save this fibre stretch
    mLambda[currentQuadPointGlobalIndex] = currentFibreStretch;

    // compute dlam/dt
    double dlam_dt = (currentFibreStretch-mLambdaLastTimeStep[currentQuadPointGlobalIndex])/(this->mNextTime-this->mCurrentTime);

    NhsSystemWithImplicitSolver& r_contraction_model = mCellMechSystems[currentQuadPointGlobalIndex];

    // Set this stretch and stretch rate
    r_contraction_model.SetStretchAndStretchRate(currentFibreStretch, dlam_dt);

    // Call RunDoNotUpdate() on the contraction model to solve it using this stretch, and get the active tension
    try
    {
        r_contraction_model.RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        rActiveTension = r_contraction_model.GetNextActiveTension();
    }
    catch (Exception& e)
    {
        #define COVERAGE_IGNORE
        // if this failed during assembling the Jacobian this is a fatal error. 
        if(assembleJacobian)
        {
            // probably shouldn't be able to get here
            EXCEPTION("Failure in solving contraction models using current stretches for assembling Jacobian");
        }
        // if this failed during assembling the residual, the stretches are too large, so we just
        // set the active tension to infinity so that the residual will be infinite
        rActiveTension = DBL_MAX;
        assert(0); // just to see if we get here, can be removed..
        return;
        #undef COVERAGE_IGNORE
    }

    // if assembling the Jacobian, numerically evaluate dTa/dlam & dTa/d(lamdot)
    if(assembleJacobian)
    {
        // get active tension for (lam+h,dlamdt)
        double h1 = std::max(1e-6, currentFibreStretch/100);
        r_contraction_model.SetStretchAndStretchRate(currentFibreStretch+h1, dlam_dt);
        r_contraction_model.RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        double active_tension_at_lam_plus_h = r_contraction_model.GetNextActiveTension();

        // get active tension for (lam,dlamdt+h)
        double h2 = std::max(1e-6, dlam_dt/100);
        r_contraction_model.SetStretchAndStretchRate(currentFibreStretch, dlam_dt+h2);
        r_contraction_model.RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        double active_tension_at_dlamdt_plus_h = r_contraction_model.GetNextActiveTension();

        rDerivActiveTensionWrtLambda = (active_tension_at_lam_plus_h - rActiveTension)/h1;
        rDerivActiveTensionWrtDLambdaDt = (active_tension_at_dlamdt_plus_h - rActiveTension)/h2;
    }

    // Re-set the stretch and stretch rate and recompute the active tension so that
    // if this guess turns out to the solution, we can just update the state variables
    //  -- not needed as AssembleSystem(true,false) [ie assemble residual] is
    //     called in ImplicitCardiacMechanicsAssembler<DIM>::Solve() above
    //     after the solve and before the update.
    //  -- The SetActiveTensionInitialGuess() would make this very fast
    //     (compared to AssembleSystem(true,false) above), but the NHS class uses the last
    //     active tension as the initial guess anyway..
    //r_contraction_model.SetStretchAndStretchRate(currentFibreStretch, dlam_dt);
    //r_contraction_model.SetActiveTensionInitialGuess(rActiveTension);
    //r_contraction_model.RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
    //assert( fabs(r_contraction_model.GetNextActiveTension()-rActiveTension)<1e-8);
}    



template class ImplicitCardiacMechanicsAssembler<2>;
template class ImplicitCardiacMechanicsAssembler<3>;


