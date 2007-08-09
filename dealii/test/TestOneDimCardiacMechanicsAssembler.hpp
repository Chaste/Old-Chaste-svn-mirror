#ifndef TEST1DCARDIACMECHANICSASSEMBLER_HPP_
#define TEST1DCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "OneDimCardiacMechanicsAssembler.hpp"

class TestOneDimCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void Testsdffd()
    {
        Triangulation<1> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(7);
        
        OneDimCardiacMechanicsAssembler mechanics(&mesh);
                                                       
        mechanics.Solve();
    }
};
#endif /*TEST1DCARDIACMECHANICSASSEMBLER_HPP_*/
