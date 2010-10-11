/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTWNTCONCENTRATION_HPP_
#define TESTWNTCONCENTRATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "MeshBasedCellPopulation.hpp"
#include "WntCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestWntConcentration : public AbstractCellBasedTestSuite
{
public:

    void TestNoWnt() throw(Exception)
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);
        p_wnt->SetType(NONE);
        p_wnt->SetCryptLength(22.0);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);   // NONE does not register as a set up Wnt Gradient (so stem cells are not moved)

        TS_ASSERT_EQUALS(p_wnt->GetType(), NONE);

        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        c_vector<double,2> location;
        location[0] = 1.5;
        location[1] = 2.3;

        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[1], 0.0, 1e-12);

        WntConcentration<2>::Destroy();
    }


    void TestLinearWntConcentration() throw(Exception)
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);
        p_wnt->SetType(LINEAR);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);

        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 1.0-height/crypt_length, 1e-9);

        // Test GetWntGradient() method
        c_vector<double,2> location;
        location[0] = 1.5;
        location[1] = 2.3;

        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[0], 0.0, 1e-12);
        // This should be equal to -1/22 = -0.0454
        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[1], -0.0454, 1e-4);

        WntConcentration<2>::Destroy();
    }

    void TestExponentialWntConcentration() throw(Exception)
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);
        p_wnt->SetType(EXPONENTIAL);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);

        TS_ASSERT_DELTA(p_wnt->GetWntConcentrationParameter(),1.0, 1e-9);

        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        // For heights above the top of the crypt (no Wnt)
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;    // for cells very close to 0 on negative side (very strong Wnt = 1.0);
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        // For normal 'in range' Wnt height
        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, exp(-height/crypt_length), 1e-9);

        // For a change in lambda
        p_wnt->SetWntConcentrationParameter(0.5);
        TS_ASSERT_DELTA(p_wnt->GetWntConcentrationParameter(),0.5, 1e-9);
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, exp(-(height/crypt_length)/0.5), 1e-9);

        // Test GetWntGradient() method

        c_vector<double,2> location;
        location[0] = 1.5;
        location[1] = 2.3;

        TS_ASSERT_THROWS_THIS(p_wnt->GetWntGradient(location)[0],
                              "No method to calculate gradient of this Wnt type");

        WntConcentration<2>::Destroy();
    }


    void TestOffsetLinearWntConcentration() throw(Exception)
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(LINEAR);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);
        p_wnt->SetWntConcentrationParameter(1.0/3.0);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);

        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        // under a third of the way up the crypt.
        height = 7.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0 - height/((1.0/3.0)*crypt_length), 1e-9);

        // more than a third of the way up the crypt.
        height = 10.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        WntConcentration<2>::Destroy();
    }

    void TestRadialWntConcentration() throw(Exception)
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(RADIAL);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);   // only fully set up when a cell population is assigned

        // Test GetWntLevel(double) method
        double height = 100;
        double wnt_level = 0.0;

        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-4);

        height = -1e-12;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-4);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0454, 1e-4);

        height = 7.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.6818, 1e-4);

        // Test GetWntLevel(CellPtr) method

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Translate mesh so that its centre is at (0,0)
        mesh.Translate(-0.5,-0.5);

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create a crypt
        MeshBasedCellPopulation<2> crypt(mesh, cells);
        p_wnt->SetCellPopulation(crypt);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);    // fully set up now

        WntConcentration<2>::Destroy();

        WntConcentration<2>::Instance()->SetType(NONE);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        crypt_length = 1.0;
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), false);    // not fully set up now it is a NONE type

        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);

        p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetCryptLength(crypt_length);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);    // set up again

        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();

        double wnt_at_cell0 = p_wnt->GetWntLevel(*cell_iter);

        double a = p_wnt->GetCryptProjectionParameterA();
        double b = p_wnt->GetCryptProjectionParameterB();
        TS_ASSERT_DELTA(a, 0.5, 1e-12);
		TS_ASSERT_DELTA(b, 2.0, 1e-12);

        while (cell_iter != crypt.End())
        {
            TS_ASSERT_DELTA(p_wnt->GetWntLevel(*cell_iter), wnt_at_cell0, 1e-12);

            // Test GetWntGradient(CellPtr) method
            c_vector<double,2> cell_location = crypt.GetLocationOfCellCentre(*cell_iter);
            double r = norm_2(cell_location);

            c_vector<double,2> expected_wnt_gradient;
            expected_wnt_gradient[0] = -cell_location[0]*pow(r,b-1.0)/(a*r);
            expected_wnt_gradient[1] = -cell_location[1]*pow(r,b-1.0)/(a*r);

            TS_ASSERT_DELTA(p_wnt->GetWntGradient(*cell_iter)[0],expected_wnt_gradient[0],1e-6);
            TS_ASSERT_DELTA(p_wnt->GetWntGradient(*cell_iter)[1],expected_wnt_gradient[1],1e-6);

            ++cell_iter;
        }

        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        WntConcentration<2>::Destroy();
    }


    void TestArchiveWntConcentration()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_grad.arch";

        // Create an output archive
        {
            WntConcentration<2>* p_wnt1 = WntConcentration<2>::Instance();
            p_wnt1->SetType(LINEAR);
            double crypt_length = 22.0;
            p_wnt1->SetCryptLength(crypt_length);
            p_wnt1->SetCryptProjectionParameterA(3.3);
            p_wnt1->SetCryptProjectionParameterB(4.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << static_cast<const WntConcentration<2>&> (*WntConcentration<2>::Instance());

            WntConcentration<2>::Destroy();
        }

        {
            WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> *p_wnt;

            double height = 21.0;
            double wnt_level = p_wnt->GetWntLevel(height);

            TS_ASSERT_DELTA(wnt_level, 1.0-height/p_wnt->GetCryptLength(), 1e-9);
            TS_ASSERT_DELTA(p_wnt->GetCryptLength(), 22.0, 1e-12);
            TS_ASSERT_DELTA(p_wnt->GetCryptProjectionParameterA(), 3.3, 1e-12);
            TS_ASSERT_DELTA(p_wnt->GetCryptProjectionParameterB(), 4.4, 1e-12);
        }

        WntConcentration<2>::Destroy();
    }


    void TestSingletonnessOfWntConcentration()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(NONE);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);

        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        TS_ASSERT_THROWS_THIS(p_wnt->SetType(NONE),"Destroy has not been called");
        TS_ASSERT_THROWS_THIS(p_wnt->SetCryptLength(10.0),"Destroy has not been called");
        WntConcentration<2>::Destroy();

        p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(LINEAR);
        p_wnt->SetCryptLength(crypt_length);

        height = 100;
        wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 1.0-height/crypt_length, 1e-9);

        TS_ASSERT_THROWS_THIS(p_wnt->SetConstantWntValueForTesting(-10),"WntConcentration<DIM>::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");

        WntConcentration<2>::Destroy();
    }


    void TestWntInitialisationSetup() throw(Exception)
    {
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<WntCellCycleModel*> models;

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
            models.push_back(p_model);
        }

        // Create the crypt
        MeshBasedCellPopulation<2> crypt(mesh, cells);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        double crypt_length = 1.0;
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);

        // As there is no cell-based simulation we must explicitly initialise the cells
        crypt.InitialiseCells();

        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
        	WntCellCycleModel* p_model = static_cast<WntCellCycleModel*>(cell_iter->GetCellCycleModel());
            std::vector<double> proteins = p_model->GetProteinConcentrations();

            if (crypt.GetLocationOfCellCentre(*cell_iter)[1]==0.0)
            {
                TS_ASSERT_DELTA(proteins[5], 4.975124378109454e-03, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7], 6.002649406788524e-01, 1e-3);
                TS_ASSERT_DELTA(proteins[8], 1.00, 1e-3);
            }
            else
            {
                TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 1.0, 1e-12);
                TS_ASSERT_DELTA(proteins[5], 1.000, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7], 0.0074, 1e-3);
                TS_ASSERT_DELTA(proteins[8], 0.00, 1e-3);
            }
        }

        WntConcentration<2>::Destroy();
    }


    void TestCryptProjectionParameterAAndBGettersAndSetters()
    {
        WntConcentration<2>* p_wnt1 = WntConcentration<2>::Instance();

        p_wnt1->SetCryptProjectionParameterA(0.8);
        p_wnt1->SetCryptProjectionParameterB(1.3);

        WntConcentration<2>* p_wnt2 = WntConcentration<2>::Instance();

        TS_ASSERT_DELTA(p_wnt2->GetCryptProjectionParameterA(), 0.8, 1e-12);
        TS_ASSERT_DELTA(p_wnt2->GetCryptProjectionParameterB(), 1.3, 1e-12);
    }

};

#endif /*TESTWNTCONCENTRATION_HPP_*/
