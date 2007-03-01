#ifndef DYNAMICFINITEELASTICITYASSEMBLER_HPP_
#define DYNAMICFINITEELASTICITYASSEMBLER_HPP_


// speed notes:
//  - don't use Tensor<4,DIM>, use double[][][][] - factor of about 2 difference
//  - compile using OPTIMISATION!!!!!! - factor of about 200 difference when 
//     assembling system!!
//  - make reused variables static rather then repeatedly creating them? or member 
//     variables. static may cause problems if a 2d sim is run then a 3d sim is run? 
//     or probably not.


// todos: see finite elas assembler

#include "FiniteElasticityAssembler.cpp"   //cpp!

template<int DIM>
class DynamicFiniteElasticityAssembler : public FiniteElasticityAssembler<DIM>
{
private:

/* Inherited from base class:
    Triangulation<DIM>*  mpMesh;
    FESystem<DIM>        mFeSystem;
    DoFHandler<DIM>      mDofHandler;

    SparsityPattern      mSparsityPattern;

    SparseMatrix<double> mJacobianMatrix;
    Vector<double>       mCurrentSolution;
    Vector<double>       mResidual;
    
    std::string          mOutputDirectoryFullPath;

    double dTdE[DIM][DIM][DIM][DIM]; 

    AbstractIncompressibleMaterialLaw<DIM>*  mpMaterialLaw;
    Vector<double>       mBodyForce;
    double               mDensity;
    const unsigned       PRESSURE_COMPONENT_INDEX; 

    std::map<unsigned,double> mBoundaryValues;

    void AssembleSystem(bool assembleResidual, bool assembleJacobian);
    void ApplyDirichletBoundaryConditions(bool assembleResidual, bool assembleJacobian);

    double CalculateResidualNorm();
    void OutputResults(unsigned newtonIteration);
    SparseMatrix<double> mNumericalJacobianMatrix;
 */


    Vector<double> mSolutionAtLastTimestep;

    double mTstart;
    double mTend;
    double mDt;
    double mDtInverse;
    bool mTimesSet;
    
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
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
    DynamicFiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                                     AbstractIncompressibleMaterialLaw<DIM>*  pMaterialLaw,
                                     Vector<double> bodyForce,
                                     double density,
                                     std::string outputDirectory,
                                     unsigned degreeOfBasesForPosition=2,
                                     unsigned degreeOfBasesForPressure=1);
    virtual ~DynamicFiniteElasticityAssembler();
    
    /**
     *  Set the start and end times, and dt, for the simulation. Must be called before Solve()
     */
    void SetTimes(double Tstart, double Tend, double dt);

    void Solve();


/* Inherited
    void SetBoundaryValues(std::map<unsigned, double> boundary_values);
    Vector<double>& GetSolutionVector();
    DoFHandler<DIM>& GetDofHandler();
    void ComputeNumericalJacobian();
    void CompareJacobians();
*/    
};




#endif /*DYNAMICFINITEELASTICITYASSEMBLER_HPP_*/
