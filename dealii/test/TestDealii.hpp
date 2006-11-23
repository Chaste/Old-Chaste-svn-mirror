#ifndef TESTDEALII_HPP_
#define TESTDEALII_HPP_

#include <cxxtest/TestSuite.h>
#include <grid/tria.h>
//#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
//#include <grid/tria_accessor.h>
//#include <grid/tria_iterator.h>


class TestDealii : public CxxTest::TestSuite
{
public:
    void testDealii()
    {
        Triangulation<2>     triangulation;
        GridGenerator::hyper_cube (triangulation, -1, 1);
        std::cout << "hello deal ii\n";
    }
};

#endif /*TESTDEALII_HPP_*/
