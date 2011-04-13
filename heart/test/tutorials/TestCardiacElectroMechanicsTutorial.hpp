
/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_
#define TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_

/*
 * = Cardiac Electro-mechanical Problems =
 * 
 * EMPTYLINE

 * == Introduction ==
 * 
 * EMPTYLINE
 * 
 * The tutorial explains how electro-mechanics problems can be solved in Chaste. The reader should certainly read
 * the electro-physiological tutorials before this tutorial, and it is helpful to have also had a look at
 * the tutorial on solving general solid mechanics problems.
 *
 * The equations of cardiac electro-mechanics are written down in Section 4.2 of the PDF on equations and
 * finite element implementations in ChasteGuides -> Miscellaneous information. '''Note:''' By default we do
 * not these full equations: the mechanics information is not coupled back to electrics, ie by default
 * the conductivities do not depend  on deformation, and cell models do not get affected by stretch.
 * This has to be switched on if required - see comments on mechano-electric feedback below.
 *
 * EMPTYLINE
 *
 * Before going to the code, we list the sub-models/parameters that need to be set, or can be varied,
 * in electro-mechanical problems. The last four of the below are mechanics-specific.
 *  * The geometry (see note 1 below)
 *  * The region electrically stimulated
 *  * The cell model
 *  * Electro-physiological parameters (conductivity, capacitance, surface-area-to-volume ratio)
 *  * Electro-physiological timesteps: ode and pde (but not printing timestep) (see note 2 below)
 *  * Fibre directions (and maybe sheet/normal directions) (see note 3 below)
 *  * The part of the boundary that is fixed in space
 *  * The contraction model [the model which takes in electrical variables (voltage or calcium typically), and
 *  returns cellular active tension]
 *  * The material law [the strain-energy function]
 *  * Mechanics timesteps: mechanics update timestep, contraction model ode timestep. (see note 4 below)
 *
 * Notes:
 *  * ''Meshes:'' Two meshes for the geometry are required, one for the electrics solve and one for the mechanics.
 * The mechanics mesh would ideally be coarser but any two meshes are technically possible. The meshes should
 * ideally both cover exactly the same geometry (ie either mesh being contained in the other), but the meshes
 * not completely overlapping is allowed - some extrapolation of quantities will then occur.
 *  * ''The electro-physiology printing timestep:'' This is not used in electro-mechanics problems; output is
 * instead written after every mechanics solve, so effectively the mechanics update timestep is equal to
 * the printing timestep.
 *  * ''Fibres:'' In electro-physiological simulations the fibre direction is in the X-direction
 * by default, but if isotropic conductivities are used the fibre direction won't be used. In mechanics
 * solves, the fibres will always be used as it determines the direction of contraction. Sheet/normal directions
 * may be used in the material law.
 *  * ''Timesteps:'' The should-divide rules are: (a) ode_timestep should-divide pde_timestep should-divide
 *  mechanics_update_timestep and (b) contraction_model_ode_timestep should-divide mechanics_update_timestep.
 * 
 * EMPTYLINE
 * 
 * The basic includes are */
#include <cxxtest/TestSuite.h>
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudy1991.hpp"
/* The includes for the second test are */
#include "CardiacElectroMechanicsProblem.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
/*
 * EMPTYLINE
 *
 * == IMPORTANT: using HYPRE ==
 *
 * Mechanics solves involve solving a nonlinear system, which is broken down into a sequence of linear solves.
 * When running problems '''in 3D, or with more elements than in the first test below''', it is vital to change the linear
 * solver to use HYPRE, an algebraic multigrid solver. Without HYRPE, the linear solve (i) may become very very slow; or
 * (ii) may not converge, in which case the nonlinear solve will (probably) not converge. HYPRE is (currently) not a
 * pre-requisite for installing Chaste, hence this is not (currently) the default linear solver for mechanics problems,
 * although this will change in the future. HYPRE should be considered a pre-requisite for large mechanics problems.
 * You can run the first test below without HYPRE, but it is certainly recommended for the second test.
 *
 * EMPTYLINE
 *
 * To use HYRPE in mechanics solves, you need to have Petsc installed with HYPRE. However, if you followed installation
 * instructions for Chaste 2.1 or later, you probably do already have Petsc installed with HYPRE.
 *
 * EMPTYLINE
 *
 * To switch on HYPRE, open the file `pde/src/solver/AbstractNonlinearElasticitySolver` and uncomment the line
 * #define MECH_USE_HYPRE
 * near the top of the file (currently: line 53).
 *
 * EMPTYLINE
 *
 * Mechanics solves being nonlinear are expensive, so it is recommended you also use `build=GccOpt_ndebug` (when running scons)
 * on larger problems.
 *
 * EMPTYLINE
 *
 * Note: Petsc unfortunately doesn't quit if you try to use HYPRE without it being installed, but it spew lots of error messages.
 *
 */

    /* EMPTYLINE
     * 
     * == Simple 2d test ==
     * 
     * EMPTYLINE
     * 
     * This test shows how to use the `CardiacElectroMechProbRegularGeom` class, which
     * inherits from the more general class `CardiacElectroMechanicsProblem` class but
     * sets up a square or cubic geometry for you.
     */
class TestCardiacElectroMechanicsTutorial : public CxxTest::TestSuite
{
public:
    void TestCardiacElectroMechanicsExample() throw(Exception)
    {
        /* All electro-mechanics problems require a cell factory as normal. This particular
         * factory stimulates the LHS side (X=0) surface. */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        /* The `CardiacElectroMechProbRegularGeom<2>` defines an electro-mechanics problem on a square. Two meshes are created
         * internally. */
        CardiacElectroMechProbRegularGeom<2> problem(KERCHOFFS2003,  // The contraction model (see below)
                                                     0.1,  // width of square (cm) 
                                                     5,    // Number mechanics elements in each direction 
                                                     10,   // Number electrics elements in each direction  
                                                     &cell_factory,
                                                     40.0, // end time
                                                     0.01, // electrics timestep (ms) 
                                                     100,  // The number of electrics dt per mechanics update - so here mech_update=1ms 
                                                     0.01, // contraction model ode timestep 
                                                     "TestCardiacElectroMechanicsExample" /* output directory */);
        /* The contraction model chosen above is KERCHOFFS2003 (Kerchoffs, Journal of Engineering Mathematics, 2003). Other possibilities
         * are 'NHS' (Niederer, Hunter, Smith, 2006), and 'NASH2004' (Nash, Progress in biophysics and molecular biology, 2004).
         *
         * EMPTYLINE
         *
         * Two meshes are created, one with five elements in each direction for the mechanics (so 5*5*2 triangles in total),
         * and a finer one for the electrics.
         * 
         * EMPTYLINE 
         * 
         * This leaves the material law, fibres direction and fixed nodes from the list above: the material
         * law is hard-coded to the Pole-zero material law, the fibre direction is by default the X-direction, and the fixed
         * nodes are automatically set (when `CardiacElectroMechProbRegularGeom` is used) to be those satistying X=0, ie
         * the left-hand edge. We discuss how to change these in the second test.
         *
         * EMPTYLINE
         *
         * Then all we have to do is call Solve.
         */
        problem.Solve();
        
        /* Go to the output directory. There should be log file (which can be used to watch progress), and a 
         * directory for the electrics output and the mechanics output. The electrics directory is not the same
         * as when running an electrics solve: the basic HDF5 data is there but there is no there is no meshalyzer
         * output, and there is always cmgui output of the ''electrics solution downsampled onto the mechanics mesh''.
         * The deformation output directory contains the deformed solution each timestep in several simple
         * matlab-readable files, and a cmgui output directory. The latter has a script for automatically loading
         * all the results.
         * 
         * EMPTYLINE
         * 
         * Visualise the results by calling `cmgui LoadSolutions.com` in the directory
         * `TestCardiacElectroMechanicsExample/deformation/cmgui` . The electrics data can be visualised on the
         * deforming mesh by using the scene (and spectrum) editor. (See cmgui website for information on how
         * to use cmgui, but very briefy: graphics -> scene editor -> select surfaces -> add, then check 'Data'. Then
         * graphics -> Spectrum editor -> min=-90, max=50.).
         *
         * EMPTYLINE
         *
         * To observe the tissue relaxing you can re-run the simulation with an end time of more than 350ms.
         */
    }



    /* EMPTYLINE
     * 
     * == Twisting cube: 3d example with varying fibre directions ==
     * 
     * EMPTYLINE
     * 
     * The second test is a longer running 3d test - the 'dont' in the name of the test
     * means it isn't run automatically. To run, remove the 'dont'. It is worth running
     * with `build=GccOpt_ndebug`, and '''see the comments about HYPRE above.'''
     * 
     * EMPTYLINE
     * 
     * This test shows how to do 3d simulations (trivial changes), and how to use 
     * `CardiacElectroMechanicsProblem`, which requires meshes and fixed nodes to be passed 
     * in, and also how to pass in fibre directions for the mechanics mesh.
     */
    void dontTestTwistingCube() throw(Exception)
    {
        /* Cell factory as normal */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-1000*1000);

        /* Set up two meshes of 1mm by 1mm by 1mm, one a `TetrahedralMesh`
         * for the electrics solve, one a (coarser) `QuadraticMesh` for the mechanics
         * solve. */
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, 0.1/*length*/, 0.1/*width*/, 0.1/*depth*/);
    
        QuadraticMesh<3> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1, 0.1 /*as above with a different stepsize*/);

        /* We choose to fix the nodes on Z=0. For this the `NonlinearElasticityTools` class
         * is helpful. The static method called below returns all nodes for which the Z value
         * (indicated by the '2' ('0' for X, '1' for Y)) is equal to 0.0. */
        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 2, 0.0);

        /* Create the problem object, which has the same interface as the the child class used
         * in the first test, except it takes in meshes and fixed nodes (as std vectors) */
        CardiacElectroMechanicsProblem<3> problem(KERCHOFFS2003,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  50,   // end time 
                                                  0.01, // electrics timestep (ms) 
                                                  100,  // The number of electrics dt per mechanics update - so here mech_update=1ms 
                                                  1.0,  // contraction model ode timestep 
                                                  "TestCardiacElectroMech3dTwistingCube" /* output directory */);

        /* The default fibre direction is the X-direction (and the default sheet plane is the XY plane). Here we show
         * how this can be changed.
         *
         * EMPTYLINE
         *
         * Fibre files should be .ortho type files (not .axi), since the sheet direction is used in the default material
         * law (see file formats documentation if you haven't come across these files, basically .axi files specify the
         * fibre directions; .ortho the fibre sheet and normal directions). For mechanics problems, the .ortho file
         * can be used to either define the fibre information PER-ELEMENT or PER-QUADRATURE-POINT (ie all the quad points
         * in all the elements). The latter provides a higher resolution description of fibres. Here we use the latter, just
         * because it is the harder case. Tthe `true` below the problem class tells the class the fibres are defined per quad
         * point. To see how this data file was generated, see below. */
        problem.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/5by5by5_fibres_by_quadpt.orthoquad", true);

        /* `SetNoElectricsOutput` is a method that is sometimes useful with a fine electrics mesh (although in this
         * case we don't call it). */ 
        bool no_electrics = false;
        if(no_electrics)
        {
            problem.SetNoElectricsOutput();
        }

        /* Now call `Solve`. This will take a while to run, so watch progress using the log file to estimate when
         * it will finish. `build=GccOpt_ndebug` will speed this up by a factor of about 5.
         */
        problem.Solve();


        /* The way the fibre file was created is given here. After defining the mechanics mesh, do: */
        //GaussianQuadratureRule<3> quad_rule(3);
        //QuadraturePointsGroup<3> quad_points(mechanics_mesh, quad_rule);
        //std::cout << quad_points.Size() << "\n";
        //for(unsigned i=0; i<quad_points.Size(); i++)
        //{
        //    ////std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << quad_points.Get(i)(2) << " ";
        //    double x = quad_points.Get(i)(0);
        //    double theta = M_PI/3 - 10*x*2*M_PI/3; // 60 degrees when x=0, -60 when x=0.1;
        //    std::cout <<  "0 " << cos(theta)  << " " << sin(theta)
        //              << " 0 " << -sin(theta) << " " << cos(theta)
        //              << " 1 0 0\n";
        //}
        /* For creating a fibre file with fibres for each element instead, we could have done */
        //for(unsigned i=0; i<mechanics_mesh.GetNumElements(); i++)
        //{
        //    double X = mechanics_mesh.GetElement(i)->CalculateCentroid()(0);
        //    //etc
        //}

        /* The one thing we haven't shown how to change is the material law. Unfortunately this is currently
         * hardcoded (ie there is no interface to change it) to the pole-zero material law. It can be manually changed
         * by altering the file `heart/src/solver/mechanics/AbstractCardiacMechanicsSolver` - search for
         * `NashHunterPoleZeroLaw`. This issue will be fixed in the near future. */
    }
};

    /* == Mechano-electric feedback ==
     *
     * As mentioned above, by default feedback of the mechanics to the electrics is not switched on, so ''by default''
     * the conductivities will not be affected by the deformation, and the stretch is not passed back to the cell-models
     * to allow for stretch-activated channels (SAC). To allow for these two features, call
     */
    //problem.UseMechanoElectricFeedback();
    /* before calling `problem.Solve()`. Note that (i) the electrics solve will slow down, since the linear system matrix now
     * varies with time (as conductivities depend on deformation), and has to be recomputed after every mechanics update; and
     * (ii) if you want a cell model that includes SAC you have to implement one. There is a single example of this in
     * the code base at the moment, see  `heart/src/odes/ionicmodels/NobleVargheseKohlNoble1998WithSac.hpp`.
     *
     * EMPTYLINE
     *
     * Further functionality and examples using M.E.F. will be added in the near future.
     *
     * EMPTYLINE
     *
     * == Other comments ==
     *
     * If you would like to apply a traction boundary condition, see the solid mechanics tutorial on how to apply tractions given
     * a `NonlinearElasticitySolver`, and then note that you can access this solver in the tests above by doing, for example:
     */
    //problem.Initialise();
    //problem.GetCardiacMechanicsSolver()->SetSurfaceTractionBoundaryConditions(boundary_elems, tractions);
    /* and then calling `problem.Solve()`. */

#endif /*TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_*/

