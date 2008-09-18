/*

Copyright (C) University of Oxford, 2008

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

#ifndef ABSTRACTNONLINEARELASTICITYASSEMBLER_HPP_
#define ABSTRACTNONLINEARELASTICITYASSEMBLER_HPP_

#include <petsc.h>
#include <vector>
#include <cmath>
#include "PetscTools.hpp"
#include "LinearSystem.hpp"
#include "AbstractIncompressibleMaterialLaw2.hpp"
#include "OutputFileHandler.hpp"

template<unsigned DIM>
class AbstractNonlinearElasticityAssembler
{
protected:
    /*< Maximum absolute tolerance for newton solve  */
    static const double MAX_NEWTON_ABS_TOL = 1e-8;
    /*< Minimum absolute tolerance for newton solve  */
    static const double MIN_NEWTON_ABS_TOL = 1e-12;
    /*< Relative tolerance for newton solve  */
    static const double NEWTON_REL_TOL = 1e-4;

    /** 
     * Number of degrees of freedom (eg equal to DIM*N + M if quadratic-linear
     * bases are used, where there are N total nodes and M vertices).
     */ 
    unsigned mNumDofs;

    /**
     *  The material laws for each element. This will either be of size
     *  1 (same material law for all elements, ie homogeneous), or size
     *  num_elem.
     */
    std::vector<AbstractIncompressibleMaterialLaw2<DIM>*> mMaterialLaws;
    /**
     *  The linear system where we store all residual vectors which are calculated
     *  and the Jacobian. Note we don't actually call Solve but solve using Petsc
     *  methods explicitly (in order to easily set num restarts etc). In the future
     *  it'll be solved using the UMFPACK direct method */ 
    LinearSystem* mpLinearSystem;
    /*< Body force vector */
    c_vector<double,DIM> mBodyForce;
    /*< Mass density of the undeformed body (equal to the density of deformed body) */
    double mDensity;

    /*< Where to write output, relative to CHASTE_TESTOUTPUT */
    std::string mOutputDirectory;
    /*< All nodes (including non-vertices) which are fixed */
    std::vector<unsigned> mFixedNodes;
    /*< The displacements of those nodes with displacement boundary conditions */
    std::vector<c_vector<double,DIM> > mFixedNodeDisplacements;

    /*< Whether to write any output */
    bool mWriteOutput;
    
    /** 
     *  The current solution, in the form (in 2d)
     *  [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *  where there are N total nodes and M vertices
     */
    std::vector<double> mCurrentSolution;
    
    /** 
     *  Storage space for a 4th order tensor used in assembling the 
     *  Jacobian (to avoid repeated memory allocation)
     */
    FourthOrderTensor2<DIM> dTdE;
    
    /*< Number of newton iterations taken in last solve */
    unsigned mNumNewtonIterations;
    
    
    /*< Deformed position: mDeformedPosition[i](j) = x_j for node i */
    std::vector<c_vector<double,DIM> > mDeformedPosition;

    /** 
     *  The solution pressures. mPressures[i] = pressure at node i (ie
     *  vertex i).
     */
    std::vector<double> mPressures;
    /**
     *  The surface tractions (which should really be non-zero) 
     *  for the boundary elements in mBoundaryElements
     */
    std::vector<c_vector<double,DIM> > mSurfaceTractions;

    /*< An optionally provided (pointer to a) function, giving body force as a function of undeformed position */ 
    c_vector<double,DIM> (*mpBodyForceFunction)(c_vector<double,DIM>&);
    /** An optionally provided (pointer to a) function, giving the surface traction as a function of 
      * undeformed position
      */ 
    c_vector<double,DIM> (*mpTractionBoundaryConditionFunction)(c_vector<double,DIM>&);
    /*< Whether the functional version of the body force is being used or not */
    bool mUsingBodyForceFunction;
    /*< Whether the functional version of the surface traction is being used or not */
    bool mUsingTractionBoundaryConditionFunction;



    virtual void FormInitialGuess()=0;
    virtual void AssembleSystem(bool assembleResidual, bool assembleJacobian)=0;
    virtual std::vector<c_vector<double,DIM> >& rGetDeformedPosition()=0;

    /**
     *  Apply the dirichlet boundary conditions to the linear system
     */   
    void ApplyBoundaryConditions(bool applyToMatrix)
    {
        assert(mFixedNodeDisplacements.size()==mFixedNodes.size());
        
        // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values
        // on the boundary nodes. However:
        // The boundary conditions on the LINEAR SYSTEM  Ju=f, where J is the
        // u the negative update vector and f is the residual is
        // u=current_soln-boundary_values on the boundary nodes
        for(unsigned i=0; i<mFixedNodes.size(); i++)
        {
            unsigned node_index = mFixedNodes[i];
            for(unsigned j=0; j<DIM; j++)
            {
                unsigned dof_index = DIM*node_index+j;
                double value = mCurrentSolution[dof_index] - mFixedNodeDisplacements[i](j);
                if (applyToMatrix)
                {
                    mpLinearSystem->ZeroMatrixRow(dof_index);
                    mpLinearSystem->SetMatrixElement(dof_index,dof_index,1);
                }
                mpLinearSystem->SetRhsVectorElement(dof_index, value);
            }
        }
    }

    /*< Calculate |r|_2 / length(r), where r is the current residual vector */
    double CalculateResidualNorm()
    {
        double norm;
        VecNorm(mpLinearSystem->rGetRhsVector(), NORM_2, &norm);
        return norm/mNumDofs;
    }

    /** 
     *  Take one newton step, by solving the linear system -Ju=f, (J the jacobian, f
     *  the residual, u the update), and picking s such that a_new = a_old + su (a
     *  the current solution) such |f(a)| is the smallest.
     */
    void TakeNewtonStep()
    {
        // compute Jacobian
        //Timer::Reset();
        AssembleSystem(true, true);
        //Timer::PrintAndReset("AssembleSystem");

        //Vec solution = mpLinearSystem->Solve();
        //Timer::PrintAndReset("Solve");

        // solve explicity with Petsc's GMRES method.
        KSP solver;
        Vec solution;
        VecDuplicate(mpLinearSystem->rGetRhsVector(),&solution);
        Mat& r_jac = mpLinearSystem->rGetLhsMatrix();

        KSPCreate(MPI_COMM_SELF,&solver);
        KSPSetOperators(solver, r_jac, r_jac, SAME_NONZERO_PATTERN);

        // set max iterations
        KSPSetTolerances(solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
        KSPSetType(solver,KSPGMRES);
        KSPGMRESSetRestart(solver,100); // gmres num restarts

        KSPSetFromOptions(solver);
        KSPSetUp(solver);

        KSPSolve(solver,mpLinearSystem->rGetRhsVector(),solution);

        ReplicatableVector update(solution);

        std::vector<double> old_solution(mNumDofs);
        for(unsigned i=0; i<mNumDofs; i++)
        {
            old_solution[i] = mCurrentSolution[i];
        }

        double best_norm_resid = DBL_MAX;
        double best_damping_value = 0.0;

        std::vector<double> damping_values;
        damping_values.reserve(12);
        damping_values.push_back(0.0);
        damping_values.push_back(0.05);
        for (unsigned i=1; i<=10; i++)
        {
            damping_values.push_back((double)i/10.0);
        }

        for (unsigned i=0; i<damping_values.size(); i++)
        {
            for(unsigned j=0; j<mNumDofs; j++)
            {
                mCurrentSolution[j] = old_solution[j] - damping_values[i]*update[j];
            }

            // compute residual
            AssembleSystem(true, false);
            double norm_resid = CalculateResidualNorm();

            std::cout << "\tTesting s = " << damping_values[i] << ", |f| = " << norm_resid << "\n" << std::flush;
            if (norm_resid < best_norm_resid)
            {
                best_norm_resid = norm_resid;
                best_damping_value = damping_values[i];
            }
        }

        if (best_damping_value == 0.0)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Residual does not decrease in newton direction, quitting");
            #undef COVERAGE_IGNORE
        }
        else
        {
            std::cout << "\tBest s = " << best_damping_value << "\n"  << std::flush;
        }

        // implement best update and recalculate residual
        for(unsigned j=0; j<mNumDofs; j++)
        {
            mCurrentSolution[j] = old_solution[j] - best_damping_value*update[j];
        }
        
        VecDestroy(solution);
        KSPDestroy(solver);
    }


public:
    AbstractNonlinearElasticityAssembler(unsigned numDofs,
                                         AbstractIncompressibleMaterialLaw2<DIM>* pMaterialLaw,
                                         c_vector<double,DIM> bodyForce,
                                         double density,
                                         std::string outputDirectory,
                                         std::vector<unsigned>& fixedNodes)
        : mNumDofs(numDofs),
          mBodyForce(bodyForce),
          mDensity(density),
          mOutputDirectory(outputDirectory),
          mFixedNodes(fixedNodes),
          mUsingBodyForceFunction(false),
          mUsingTractionBoundaryConditionFunction(false)
    {
        assert(pMaterialLaw != NULL);
        mMaterialLaws.push_back(pMaterialLaw);

        assert(DIM==2 || DIM==3);
        assert(density > 0);
        assert(fixedNodes.size()>0);
        mWriteOutput = (mOutputDirectory != "");
        
        mpLinearSystem = new LinearSystem(mNumDofs);
    }


    AbstractNonlinearElasticityAssembler(unsigned numDofs,
                                         std::vector<AbstractIncompressibleMaterialLaw2<DIM>*>& rMaterialLaws,
                                         c_vector<double,DIM> bodyForce,
                                         double density,
                                         std::string outputDirectory,
                                         std::vector<unsigned>& fixedNodes)
        : mNumDofs(numDofs),
          mBodyForce(bodyForce),
          mDensity(density),
          mOutputDirectory(outputDirectory),
          mFixedNodes(fixedNodes),
          mUsingBodyForceFunction(false),
          mUsingTractionBoundaryConditionFunction(false)
    {
        mMaterialLaws.resize(rMaterialLaws.size(), NULL);
        for(unsigned i=0; i<mMaterialLaws.size(); i++)
        {
            assert(rMaterialLaws[i] != NULL);
            mMaterialLaws[i] = rMaterialLaws[i];
        }

        assert(DIM==2 || DIM==3);
        assert(density > 0);
        assert(fixedNodes.size()>0);
        mWriteOutput = (mOutputDirectory != "");
        
        mpLinearSystem = new LinearSystem(mNumDofs);
    }


    virtual ~AbstractNonlinearElasticityAssembler()
    {
        delete mpLinearSystem;
    }
    

    /** 
     *  Solve the problem
     */    
    void Solve()
    {
        if(mWriteOutput)
        {
            WriteOutput(0);
        }
    
        // compute residual
        AssembleSystem(true, false);
        double norm_resid = this->CalculateResidualNorm();
        std::cout << "\nNorm of residual is " << norm_resid << "\n";
    
        mNumNewtonIterations = 0;
        unsigned counter = 1;
    
        // use the larger of the tolerances formed from the absolute or
        // relative possibilities
        double tol = NEWTON_REL_TOL*norm_resid;

        #define COVERAGE_IGNORE // not going to have tests in cts for everything
        if(tol > MAX_NEWTON_ABS_TOL)
        {
            tol = MAX_NEWTON_ABS_TOL;
        }
        if(tol < MIN_NEWTON_ABS_TOL)
        {
            tol = MIN_NEWTON_ABS_TOL;
        }
        #undef COVERAGE_IGNORE

        std::cout << "Solving with tolerance " << tol << "\n";
    
        while (norm_resid > tol)
        {
            std::cout <<  "\n-------------------\n"
                      <<   "Newton iteration " << counter
                      << ":\n-------------------\n";
    
            TakeNewtonStep();

            AssembleSystem(true, false);
            norm_resid = CalculateResidualNorm();
    
            std::cout << "Norm of residual is " << norm_resid << "\n";    
            if(mWriteOutput)
            {
                WriteOutput(counter);
            }
    
            mNumNewtonIterations = counter;
    
            counter++;
            if (counter==20)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Not converged after 20 newton iterations, quitting");
                #undef COVERAGE_IGNORE
            }
        }

        if (norm_resid > tol)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Failed to converge");
            #undef COVERAGE_IGNORE
        }
    
        // we have solved for a deformation so note this
        //mADeformedHasBeenSolved = true;
            
    }
    


    /**
     *  Write the current solution for the file outputdir/solution_[counter].nodes
     */
    void WriteOutput(unsigned counter)
    {
        // only write output if the flag mWriteOutput has been set
        if (!mWriteOutput)
        {
            return;
        }

        OutputFileHandler output_file_handler(mOutputDirectory, (counter==0));
        std::stringstream file_name;
        file_name << "/solution_" << counter << ".nodes";

        out_stream p_file = output_file_handler.OpenOutputFile(file_name.str());

        std::vector<c_vector<double,DIM> >& r_deformed_position = rGetDeformedPosition();
        for(unsigned i=0; i<r_deformed_position.size(); i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                *p_file << r_deformed_position[i](j) << " ";
            } 
            *p_file << "\n";
        }
        p_file->close();
    }
    
    unsigned GetNumNewtonIterations()
    {
        return mNumNewtonIterations;
    }

    /**
      * Set a function which gives body force as a function of X (undeformed position)
      * Whatever body force was provided in the constructor will now be ignored
      */
    void SetFunctionalBodyForce(c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>&))
    {
        mUsingBodyForceFunction = true;
        mpBodyForceFunction = pFunction;
    }
};

#endif /*ABSTRACTNONLINEARELASTICITYASSEMBLER_HPP_*/
