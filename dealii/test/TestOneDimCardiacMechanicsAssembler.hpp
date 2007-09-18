#ifndef TEST1DCARDIACMECHANICSASSEMBLER_HPP_
#define TEST1DCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "OneDimCardiacMechanicsAssembler.hpp"
#include "ImplicitOneDimCardiacMechanicsAssembler.hpp"

class TestOneDimCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestUncoupled() throw(Exception)
    {
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(7);
        
        OneDimCardiacMechanicsAssembler mechanics(&mesh);

        std::vector<double> active_tension(mechanics.GetTotalNumQuadPoints(), 0.5);
        
        mechanics.SetActiveTension( active_tension );
        mechanics.Solve();
        
        std::vector<Vector<double> > undeformed_position = mechanics.rGetUndeformedPosition();
        std::vector<Vector<double> > deformed_position = mechanics.rGetDeformedPosition();
        
        // not i=0 as X=0
        double factor = deformed_position[0](1)/undeformed_position[0](1);
        
        TS_ASSERT_LESS_THAN(factor, 1.0);
        
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            if(undeformed_position[0](i)!=0)
            {
                TS_ASSERT_DELTA(deformed_position[0](i)/undeformed_position[0](i), factor, 1e-4);
            }
        }
    }
    

    void TestUncoupledImplicit() throw(Exception)
    {
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(7);
        
        ImplicitOneDimCardiacMechanicsAssembler mechanics(&mesh);
        std::vector<double> caI(mechanics.GetTotalNumQuadPoints(), 0.02);
        
        mechanics.SetIntracellularCalciumConcentration( caI );
        mechanics.Solve(0, 0.1, 0.01);
        
        std::vector<Vector<double> > undeformed_position = mechanics.rGetUndeformedPosition();
        std::vector<Vector<double> > deformed_position = mechanics.rGetDeformedPosition();
        
        // not i=0 as X=0
        double factor = deformed_position[0](1)/undeformed_position[0](1);
        
        TS_ASSERT_LESS_THAN(factor, 1.0);
        
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            if(undeformed_position[0](i)!=0)
            {
                TS_ASSERT_DELTA(deformed_position[0](i)/undeformed_position[0](i), factor, 1e-4);
            }
        }
    }
    
};
#endif /*TEST1DCARDIACMECHANICSASSEMBLER_HPP_*/
