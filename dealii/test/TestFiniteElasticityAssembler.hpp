#ifndef TESTFINITEELASTICITYASSEMBLER_HPP_
#define TESTFINITEELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"


#define DIMENSION 2

class TestFiniteElasticityAssembler : public CxxTest::TestSuite
{
public:
    void testFiniteElasticityAssembler() throw(Exception)
    {
        Vector<double> body_force(DIMENSION);
        body_force(0) = 6.0;
    
        MooneyRivlinMaterialLaw<DIMENSION> mooney_rivlin_law(2.0,2.0);


        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);

        FiniteElasticity<DIMENSION> finite_elasticity(&mesh, 
                                                      &mooney_rivlin_law,
                                                      body_force,
                                                      1.0,
                                                      "finite_elas/simple");
        finite_elasticity.Solve();

//    ExponentialMaterialLaw<3> exponential_law(2.0,1.1);
//    FiniteElasticity<3> finite_elasticity2(&exponential_law, body_force, 1.0);
//    finite_elasticity2.Solve();

    }
};
#endif /*TESTFINITEELASTICITYASSEMBLER_HPP_*/
