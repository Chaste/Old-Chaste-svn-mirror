#ifndef FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_
#define FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_

#include "FiniteElasticityAssembler.cpp"
#include "AbstractGrowingTumourSourceModel.hpp"
#include "GrowthByConstantMassOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"

template<int DIM>
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
    std::vector<GrowthByConstantMassOdeSystem<DIM>*> mpGrowthOdeSystems;
    EulerIvpOdeSolver mOdeSolver;

    AbstractGrowingTumourSourceModel<DIM>* mpSourceModel;



    /**
     *  Reimplemented to include growth term (only a minor change)
     */
    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix,
                                   bool                                            assembleResidual,
                                   bool                                            assembleJacobian);

    
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
                                        unsigned degreeOfBasesForPosition=2,
                                        unsigned degreeOfBasesForPressure=1);


    virtual ~FiniteElasticityAssemblerWithGrowth();

    void SetTimes(double Tstart, double Tend, double odeDt);
    void Run();

/* Inherited    
    virtual void Solve();
    Vector<double>& GetSolutionVector();
    DoFHandler<DIM>& GetDofHandler();
    void ComputeNumericalJacobian();
    void CompareJacobians();
*/    
};




#endif /*FINITEELASTICITYASSEMBLERWITHGROWTH_HPP_*/
