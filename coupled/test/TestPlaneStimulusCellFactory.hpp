#ifndef TESTPLANESTIMULUSCELLFACTORY_HPP_
#define TESTPLANESTIMULUSCELLFACTORY_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "PlaneStimulusCellFactory.hpp"

class TestPlaneStimulusCellFactory : public CxxTest::TestSuite
{
public:
    void TestBasicContructor() throw (Exception)
    {
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);
        PlaneStimulusCellFactory<3> cell_factory1;
        PlaneStimulusCellFactory<3> cell_factory2(0.001, -100); // ode time step & stimulus voltage
        cell_factory1.SetMesh(&mesh);
        cell_factory2.SetMesh(&mesh);
        
        for (unsigned node_num=0; node_num<mesh.GetNumNodes(); node_num++)
        {
            Node<3>* p_node=mesh.GetNode(node_num);
            AbstractCardiacCell* p_cell1=cell_factory1.CreateCardiacCellForNode(node_num);
            AbstractCardiacCell* p_cell2=cell_factory2.CreateCardiacCellForNode(node_num);
            // compute 1 second
            p_cell1->Compute(0.0,1.0);
            p_cell2->Compute(0.0,1.0);
            if (p_node->GetPoint()[0] == 0.0)
            {
                // stimulated cell
                TS_ASSERT_DELTA(p_cell1->GetVoltage(), 161.1, 0.1);
                TS_ASSERT_DELTA(p_cell2->GetVoltage(), 42.0, 0.1);
            }
            else
            {
                // unstimulated cell
                TS_ASSERT_DELTA(p_cell1->GetVoltage(), -83.8, 0.1);
                TS_ASSERT_DELTA(p_cell2->GetVoltage(), -83.8, 0.1);
            }
            delete p_cell1;
            delete p_cell2;
        }
        
    }

};

#endif /*TESTPLANESTIMULUSCELLFACTORY_HPP_*/
