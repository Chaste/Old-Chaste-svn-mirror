/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_
#define FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_

#include "FiniteElasticityAssembler.cpp"
#include "AbstractGrowingTumourSourceModel.hpp"
#include "GrowthByConstantMassOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"

//todos:
// doxygen
// writebasicoutput method - move outside?

// use something like simulation time?


#define NON_GROWING_REGION 98
#define GROWING_REGION 99


/**
 *  INSERT COMMENTS
 */
template<unsigned DIM>
class FiniteElasticityAssemblerWithGrowth : public FiniteElasticityAssembler<DIM>
{
protected:
    /* Inherited:
        Triangulation<DIM>*  mpMesh;
        FESystem<DIM>        mFeSystem;
        DoFHandler<DIM>      mDofHandler;
        SparsityPattern      mSparsityPattern;
        SparseMatrix<double> mJacobianMatrix;
        Vector<double>       mCurrentSolution;
        Vector<double>       mResidual;
        std::string          mOutputDirectoryFullPath;
        bool                 mWriteOutput;
        FourthOrderTensor<DIM> dTdE;
        AbstractIncompressibleMaterialLaw<DIM>*  mpMaterialLaw;
        Vector<double>       mBodyForce;
        double               mDensity;
        const unsigned       PRESSURE_COMPONENT_INDEX;
        std::map<unsigned,double> mBoundaryValues;
        void AssembleSystem(bool assembleResidual, bool assembleJacobian);
        void ApplyDirichletBoundaryConditions(bool assembleResidual, bool assembleJacobian);
        void OutputResults(unsigned counter);
        void TakeNewtonStep();
    */
    
    double mTstart;
    double mTend;
    double mOdeDt;
    bool   mTimesSet;
    
    bool mNoRefinement;
    
    Vector<double> mGrowthValuesAtVertices;
    std::vector<GrowthByConstantMassOdeSystem<DIM>*> mGrowthOdeSystems;
    EulerIvpOdeSolver mOdeSolver;
    
    AbstractGrowingTumourSourceModel<DIM>* mpSourceModel;
    
    double mAverageElementVolume;
    
    /**
     *  Reimplemented to include growth term (only a minor change)
     */
    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix,
                                   bool                                            assembleResidual,
                                   bool                                            assembleJacobian);

    virtual void WriteStresses(unsigned counter);
    
    void WriteGrowthValuesAtVertices(unsigned counter);
    
    bool RefineOvergrownElements();
                                   
public:
    /**
     *  Constructor
     *  
     *  @param pMesh A pointer to a dealii mesh. Note, the mesh must have some surface
     *   elements which have had their boundary indicator set to FIXED_BOUNDARY
     *  @param pMaterialLaw A pointer to an incompressible material law
     *  @bodyForce A vector of size DIM represents the body force (force per unit volume)
     *  @density The mass density. Must be strictly positive
     *  @outputDirectory The output directory, relative the the testoutput directory. If 
     *   empty no output is written
     *  @degreeOfBasesForPosition Degree of the polynomials used for interpolating positions.
     *   Defaults to 2, ie quadratic interpolation 
     *  @degreeOfBasesForPressure Degree of the polynomials used for interpolating pressue.
     *   Defaults to 2, ie linear interpolation
     */
    FiniteElasticityAssemblerWithGrowth(Triangulation<DIM>* pMesh,
                                        AbstractIncompressibleMaterialLaw<DIM>*  pMaterialLaw,
                                        Vector<double> bodyForce,
                                        double density,
                                        std::string outputDirectory,
                                        AbstractGrowingTumourSourceModel<DIM>* pSourceModel,
                                        unsigned degreeOfBasesForPosition=2,
                                        unsigned degreeOfBasesForPressure=1);
                                        
                                        
    virtual ~FiniteElasticityAssemblerWithGrowth();
    
    void SetTimes(double Tstart, double Tend, double odeDt);
    void Run();
    
    
    /**
     *  Returns true if there is an growth ode system associated with
     *  this vertex in the mesh. Mainly for testing purposes
     */
    bool IsGrowingNode(unsigned vertexIndex);
    
    /** 
     *  Switches off refinement/coarsening of overgrown/shrunken elements
     */
    void SetNoRefinement();
    
    /* Inherited
        virtual void StaticSolve(bool writeOutput=true);
        DoFHandler<DIM>& rGetDofHandler();
        void ComputeNumericalJacobian();
        void CompareJacobians();
    */
};




#endif /*FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_*/
