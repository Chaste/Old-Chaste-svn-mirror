#ifndef TESTFAILINGDIVISION3D_HPP_
#define TESTFAILINGDIVISION3D_HPP_


class TestFailingDivision3d : public AbstractCellBasedTestSuite
{
public:

    /* This test will pass if you run it for 10 hours, but fails for longer times when a cell tries to divide
     *
     */
    void TestFailsWhenCellDivides() throw (Exception)
    {
        unsigned width = 4;
        unsigned height = 4;
        unsigned depth = 1;

        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = p_mesh->GetNumAllNodes();

        unsigned num_epithelial_cells = (width+1)*(height+1);

    	boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        std::vector<CellPtr> cells;
		for (unsigned i=0; i<num_nodes-num_epithelial_cells; i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(3);

            CellPtr p_differentiated_cell(new Cell(p_state, p_model));
			cells.push_back(p_differentiated_cell);
        }

        for (unsigned i=num_nodes-num_epithelial_cells; i<num_nodes; i++)
        {
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetDimension(3);

            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->InitialiseCellCycleModel();

			cells.push_back(p_epithelial_cell);
        }

        // Test Save with a MeshBasedCellPopulationWithGhostNodes
        MeshBasedCellPopulation<3> cell_population(*p_mesh, cells);

        CellBasedSimulation<3> simulator(cell_population);

        simulator.SetOutputDirectory("TestFailsWhenCellDivides");
        simulator.SetEndTime(20.0);			// Keep at 10 hours max and no cell divisions will occur. Cell divisions don't seem to work yet.

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        simulator.Solve();
    }

};

#endif /* TESTFAILINGDIVISION3D_HPP_ */
