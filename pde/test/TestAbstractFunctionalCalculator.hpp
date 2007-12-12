#ifndef TESTABSTRACTFUNCTIONALCALCULATOR_HPP_
#define TESTABSTRACTFUNCTIONALCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractFunctionalCalculator.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

// Returns 1.0 everywhere so that the total integral over the mesh of
// this integrand is just the volume of the mesh. For testing.
template<unsigned DIM>
class VolumeIntegralFunctionalCalculator : public AbstractFunctionalCalculator<DIM,DIM,1>
{
    double GetIntegrand(ChastePoint<DIM> &rX,
                        c_vector<double,1> &rU,
                        c_matrix<double,1,DIM> &rGradU)
    {
        return 1.0;
    }
};

class TestAbstractFunctionalCalculator : public CxxTest::TestSuite
{
public:
    void TestWithVolumeIntegralFunctionalCalculator()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        VolumeIntegralFunctionalCalculator<2> functional_calculator;
        
        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes(), 0.0);
        
        double result = functional_calculator.Calculate(mesh,vec);
        TS_ASSERT_DELTA(result, 1.0, 1e-6); // 'volume' of the square equals 1
        
        Vec bad_vec = PetscTools::CreateVec(mesh.GetNumNodes()+1, 0.0);
        TS_ASSERT_THROWS_ANYTHING(functional_calculator.Calculate(mesh,bad_vec));
        
        VecDestroy(vec);
        VecDestroy(bad_vec);
    }
};

#endif /*TESTABSTRACTFUNCTIONALCALCULATOR_HPP_*/
