#ifndef TESTDISCRETESYSTEMFORCECALCULATOR_HPP_
#define TESTDISCRETESYSTEMFORCECALCULATOR_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <vector>
#include "DiscreteSystemForceCalculator.hpp"
#include "TrianglesMeshReader.cpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"

class TestDiscreteSystemForceCalculator : public CxxTest::TestSuite
{
public:
    void TestCalculateFtAndFn() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        HoneycombMeshGenerator generator(7, 5, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
        Tissue<2> tissue(*p_mesh, cells);
        tissue.SetGhostNodes(ghost_node_indices);

        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
        
        DiscreteSystemForceCalculator calculator(meineke_spring_system); 
        
        std::vector<double> Ft_and_Fn = calculator.CalculateFtAndFn(8,1.0);
        TS_ASSERT_DELTA(Ft_and_Fn[0], 0.0, 1e-5);
        TS_ASSERT_DELTA(Ft_and_Fn[1], 0.0, 1e-5);
    }
};

#endif /*TESTDISCRETESYSTEMFORCECALCULATOR_HPP_*/
