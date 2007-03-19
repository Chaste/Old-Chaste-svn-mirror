#ifndef FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_
#define FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_

#include "FiniteElasticityAssembler.cpp"
#include "AbstractGrowingTumourSourceModel.hpp"
#include "GrowthByConstantMassOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"

//todos: 
// doxygen, 
// proper interpolation of g
// interpolate g onto the refined mesh
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
    double dTdE[DIM][DIM][DIM][DIM]; 
    AbstractIncompressibleMaterialLaw<DIM>*  mpMaterialLaw;
    Vector<double>       mBodyForce;
    double               mDensity;
    const unsigned       PRESSURE_COMPONENT_INDEX; 
    std::map<unsigned,double> mBoundaryValues;
    void AssembleSystem(bool assembleResidual, bool assembleJacobian);
    void ApplyDirichletBoundaryConditions(bool assembleResidual, bool assembleJacobian);
    void OutputResults(unsigned counter);
    double CalculateResidualNorm();
    void TakeNewtonStep();
*/

    double mTstart;
    double mTend;
    double mOdeDt;
    bool   mTimesSet;
    
    Vector<double> mGrowthValuesAtVertices;
    std::vector<GrowthByConstantMassOdeSystem<DIM>*> mGrowthOdeSystems;
    EulerIvpOdeSolver mOdeSolver;

    AbstractGrowingTumourSourceModel<DIM>* mpSourceModel;

    double mAverageElementVolume;
    
    void WriteBasicOutput(unsigned counter);

    /**
     *  Reimplemented to include growth term (only a minor change)
     */
    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix,
                                   bool                                            assembleResidual,
                                   bool                                            assembleJacobian);

    /**
     *  Refine the elements which have gotten too large
     */
    bool RefineOvergrownElements(unsigned);
    bool mUseRefinement;

    
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

    void DoNotUseRefinement();
    void UseRefinement();

/* Inherited    
    virtual void Solve();
    Vector<double>& rGetSolutionVector();
    DoFHandler<DIM>& rGetDofHandler();
    void ComputeNumericalJacobian();
    void CompareJacobians();
*/    
};




#endif /*FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_*/
