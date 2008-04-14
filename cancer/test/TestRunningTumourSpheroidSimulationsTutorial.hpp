/*
 * 
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 * 
 * 
 */  
#ifndef TESTTUMOURSPHEROIDTUTORIAL_HPP_
#define TESTTUMOURSPHEROIDTUTORIAL_HPP_
/* 
 * In this tutorial we show how Chaste is used to run discrete tumour 
 * spheroid simulations. These types of simulation are similar to
 * crypt simulations, in that they consist of cell cycle models and 
 * discrete mechanics laws determining how cells divide and move, but
 * these are coupled to a PDE determining the concentration of nutrients,
 * eg oxygen, throughout the domain. Also, unlike the crypt simulation,
 * the domain grows substantially as the simulation runs.
 * 
 * The main differences between this tutorial and the crypt tutorials are 
 * (i) a PDE is defined, to be used in the simulation, (ii) a non-periodic mesh 
 * is used, and (iii) the cell-cell force law is defined and explicitly
 * used (in the previous tutorials, the default cell-cell force law was used). 
 * 
 * EMPTYLINE
 * 
 * The first thing that needs to be done, when writing any Chaste test,
 * is to include the following header
 */
#include <cxxtest/TestSuite.h>
/* The following have to be included, for technical reasons....... */
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/* This header defines a helper class that is useful for generating a 
 * vector of cells */
#include "CellsGenerator.hpp"
/* These are the classes that will be used in these tests. 
 * {{{TissueSimulationWithNutrients}}} is used for tumour spheroid
 * simulations. 
 */
#include "TissueSimulationWithNutrients.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "CellwiseNutrientSinkPde.hpp"


/* Now we can define the test class, and write the test
 */
class TestRunningTumourSpheroidSimulationsTutorial : public CxxTest::TestSuite
{
public: 
    void TestSpheroidTutorial() throw(Exception)
    {        
        /* This first line can be ignored, it's a macro which just says 
         * don't run this test if in parallel. */
        EXIT_IF_PARALLEL; // defined in PetscTools.hpp

        /* The first thing to do, as before, is to set up the start time and 
         * reset the parameters. */
        SimulationTime::Instance()->SetStartTime(0.0);
        CancerParameters::Instance()->Reset();
        CancerParameters::Instance()->SetHepaOneParameters();

        /* Now we want to create a ''non-periodic'' 'honeycomb' mesh.
         * We use the honeycomb mesh generator, as before, saying 10 cells wide
         * and 10 cells high. Note that the thickness of the ghost nodes layer is
         * 0, ie no ghost nodes, and the {{{false}}} indicates not cylindrical.
         */
        HoneycombMeshGenerator generator(10, 10, 0, false);
        /* Get the mesh. Note we call {{{GetMesh()}}} rather than {{{GetCyclindricalMesh}}},
         * and that a {{{ConformingTetrahedralMesh}}} is returned. */
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

 
        /* Next, we need to create some cells. Unlike before, we don't just use
         * the {{{CellsGenerator}}} class, but do it manually, in a loop. First,
         * define the cells vector. */
        std::vector<TissueCell> cells;  
        /* then loop over the nodes... */        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /*.. then create a cell, and giving a particular cell cycle model 
             * - {{{SimpleOxygenBasedCellCycleModel}}}. The index of the node that 
             * this cell is related to also needs to be given. */ 
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel());
            cell.SetNodeIndex(i);
            
            /* Now, we define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a !HepaOne cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
             */ 
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  CancerParameters::Instance()->GetHepaOneCellG1Duration()
                                  + CancerParameters::Instance()->GetSG2MDuration() );
            /* .. then we set the birth time and push the cell back into the vector
             * of cells. */
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);            
        }

        /* Now that we have defined the cells, we can define the Tissue. This time it
         * is just a mesh-based tissue (ie not a {{{MeshBasedTissueWithGhostNodes()}}}. 
         * Again, the constructor takes in the mesh and the cells vector. */
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        
        /* Recall that in the Wnt based crypt simulation, we defined a singleton class
         * which cell-cycles used to get the wnt concentration. Here, we do the same kind
         * of thing, but using the singletom {{{CellwiseData}}} class, which stores the 
         * value of the current nutrient concentration, for each cell. We have to
         * tell the {{{CellwiseData}}} object how many nodes and variables per node there
         * are (in this case, 1 variable per node, ie the oxygen concentration), and 
         * the tissue.
         */
        CellwiseData<2>::Instance()->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        CellwiseData<2>::Instance()->SetTissue(tissue);
        /* Then we have to initialise the oxygen concentration for each node (to 1.0), by
         * calling {{{SetValue}}}. This takes in the concentration, and the node
         * which this concentration is for .*/
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellwiseData<2>::Instance()->SetValue(1.0, p_mesh->GetNode(i));
        }
        
        /* Next we instantiate an instance of the PDE class which we defined above. 
         * This will be passed into the simulator. The !CellwiseNutrientSinkPde is
         * a Pde class which inherits from !AbstractLinearEllipticPde, and represents
         * the PDE: u_xx + u_yy = k(x) u, where k(x) = 0.03 (the coefficient below) 
         * if x is in a live cell, and k(x)=0 if x is within a necrotic cell
         */
        CellwiseNutrientSinkPde<2> pde(tissue, 0.03);
        
        /* There are a several different cell-cell force laws possible, which can be
         * passed into the simulator. Here, we
         * create a {{{Meineke2001SpringSystem}}}, which uses a triangulation
         * to determine which cells are connected, and assumes a linear spring
         * between any connected cells. We can the method {{{UseCutoffPoint}}}
         * on the spring system before passing it into the simulator. This tells
         * it to return zero force if two cells are more than 3 units 
         * (=3 cell widths) away from each other. This is necessary when no ghost 
         * nodes are used.
         */
        Meineke2001SpringSystem<2> spring_system(tissue);
        spring_system.UseCutoffPoint(3);   
        
        /* 
         * The simulator object for these problems is 
         * {{{TissueSimulationWithNutrients}}}. We pass in the tissue, the
         * mechanics system, and the PDE.
         */ 
        TissueSimulationWithNutrients<2> simulator(tissue, &spring_system, &pde);
        
        /* As with {{{CryptSimulation2d}}} (which inherits from the same base class
         * as {{{TissueSimulationWithNutrients}}}), we can set the output directory
         * and end time. */
        simulator.SetOutputDirectory("SpheroidTutorial");
        simulator.SetEndTime(10.0);
                
        /* Solve. */
        simulator.Solve();

        /* Finally, call {{{Destroy()}}} on the singleton classes. The results
         * can be visualised as in the previous test. */
        SimulationTime::Destroy();
        CellwiseData<2>::Destroy();        
    }
};
#endif /*TESTTUMOURSPHEROIDTUTORIAL_HPP_*/
