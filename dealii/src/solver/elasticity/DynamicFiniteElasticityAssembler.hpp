#ifndef DYNAMICFINITEELASTICITYASSEMBLER_HPP_
#define DYNAMICFINITEELASTICITYASSEMBLER_HPP_


// speed notes:
//  - don't use Tensor<4,DIM>, use double[][][][] - factor of about 2 difference
//  - compile using OPTIMISATION!!!!!! - factor of about 200 difference when 
//     assembling system!!
//  - make reused variables static rather then repeatedly creating them? or member 
//     variables. static may cause problems if a 2d sim is run then a 3d sim is run? 
//     or probably not.


// todos: useful output, bcc object, neumann?, heterogeneity, nondim?

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
*/
    SparseMatrix<double> mNumericalJacobianMatrix;

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
    DynamicFiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                                     AbstractIncompressibleMaterialLaw<DIM>*  pMaterialLaw,
                                     Vector<double> bodyForce,
                                     double density,
                                     std::string outputDirectory,
                                     unsigned orderOfBasesForPosition=2,
                                     unsigned orderOfBasesForPressure=1);
    virtual ~DynamicFiniteElasticityAssembler();
    

    void SetTimes(double Tstart, double Tend, double dt);

    void Solve();

    void ComputeNumericalJacobian();
    void CompareJacobians();

/* Inherited
    void SetBoundaryValues(std::map<unsigned, double> boundary_values);
    Vector<double>& GetSolutionVector();
    DoFHandler<DIM>& GetDofHandler();
*/    
};




#endif /*DYNAMICFINITEELASTICITYASSEMBLER_HPP_*/
