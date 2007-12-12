#ifndef TESTABSTRACTFUNCTIONALCALCULATOR_HPP_
#define TESTABSTRACTFUNCTIONALCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractFunctionalCalculator.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"

// Returns 1.0 everywhere so that the total integral over the mesh of
// this integrand is just the volume of the mesh. For testing.
template<unsigned DIM>
class VolumeCalculator : public AbstractFunctionalCalculator<DIM,DIM,1>
{
    double GetIntegrand(ChastePoint<DIM> &rX,
                        c_vector<double,1> &rU,
                        c_matrix<double,1,DIM> &rGradU)
    {
        return 1.0;
    }
};

// Check x and u are interpolated correctly
class ExampleFunctionalOne : public AbstractFunctionalCalculator<2,2,2>
{
    double GetIntegrand(ChastePoint<2> &rX,
                        c_vector<double,2> &rU,
                        c_matrix<double,2,2> &rGradU)
    {
        return rX[0]*rU[0] + rX[1]*rU[1];
    }
};

// Check grad_u is interpolated correctly
class ExampleFunctionalTwo : public AbstractFunctionalCalculator<2,2,2>
{
    double GetIntegrand(ChastePoint<2> &rX,
                        c_vector<double,2> &rU,
                        c_matrix<double,2,2> &rGradU)
    {
        return rX[0]*rU[0] + rX[1]*rU[1] + 0.5*(rGradU(0,0)+rGradU(0,1)+rGradU(1,0)+rGradU(1,1));
    }
};

class TestAbstractFunctionalCalculator : public CxxTest::TestSuite
{
public:
    void TestWithVolumeCalculator()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        VolumeCalculator<2> volume_calculator;
        
        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes(), 0.0);
        
        double result = volume_calculator.Calculate(mesh,vec);
        TS_ASSERT_DELTA(result, 1.0, 1e-6); // 'volume' of the square equals 1
        
        Vec bad_vec = PetscTools::CreateVec(mesh.GetNumNodes()+1, 0.0);
        TS_ASSERT_THROWS_ANYTHING(volume_calculator.Calculate(mesh,bad_vec));
        
        VecDestroy(vec);
        VecDestroy(bad_vec);
    }

    void TestWithExampleFunctionals()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        // Test interpolation of x and u
        // Integrate x^2 + 2y over the unit square
        // = 4/3
        ExampleFunctionalOne calculator;
        
        Vec petsc_vec = PetscTools::CreateVec(2*mesh.GetNumNodes(), 2.0);
        DistributedVector::SetProblemSize(mesh.GetNumNodes());
        DistributedVector vec(petsc_vec);
        DistributedVector::Stripe u(vec, 0);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            Node<2>* p_node = mesh.GetNode(index.Global);
            u[index] = p_node->rGetLocation()[0];
        }
        vec.Restore();
        
        double result = calculator.Calculate(mesh, petsc_vec);
        TS_ASSERT_DELTA(result, 4.0/3.0, 1e-6);
        
        // Test interpolation of grad_u
        // Integrate x^2 + y^2 + 1 over the unit square
        // = 5/3
        ExampleFunctionalTwo other_calculator;

        DistributedVector::Stripe v(vec, 1);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            Node<2>* p_node = mesh.GetNode(index.Global);
            u[index] = p_node->rGetLocation()[0];
            v[index] = p_node->rGetLocation()[1];
        }
        vec.Restore();
        
        result = other_calculator.Calculate(mesh, petsc_vec);
        TS_ASSERT_DELTA(result, 1 + 2.0/3.0, 1e-6);
        
        VecDestroy(petsc_vec);
    }

};

#endif /*TESTABSTRACTFUNCTIONALCALCULATOR_HPP_*/
