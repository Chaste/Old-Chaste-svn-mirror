/*

Copyright (C) University of Oxford, 2008

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
#ifndef TESTCELLSGENERATOR_HPP_
#define TESTCELLSGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "FixedCellCycleModelCellsGenerator.hpp"
#include "SimpleWntCellCycleModelCellsGenerator.hpp"
#include "StochasticCellCycleModelCellsGenerator.hpp"
#include "StochasticWntCellCycleModelCellsGenerator.hpp"
#include "TysonNovakCellCycleModelCellsGenerator.hpp"
#include "WntCellCycleModelCellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestCellsGenerator : public AbstractCancerTestSuite
{
public:

    void TestFixedCellCycleModelCellsGeneratorGenerateBasic() throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;

        FixedCellCycleModelCellsGenerator<2> generator;
        generator.GenerateBasic(cells, mesh);

        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), -(double)(i), 1e-9);
        }
    }

    void TestFixedCellCycleModelCellsGeneratorGenerateForCrypt() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;

        std::vector<TissueCell> cells;

        double y0 = 0.2;
        double y1 = 1.0;
        double y2 = 2.0;
        double y3 = 3.0;
        
        FixedCellCycleModelCellsGenerator<2> generator;
        generator.GenerateForCrypt(cells, *p_mesh, true, y0, y1, y2 ,y3 );

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = cells[i].GetCellCycleModel()->GetGeneration();
            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }
        }
    }

    void TestStochasticCellCycleModelCellsGenerator() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;

        StochasticCellCycleModelCellsGenerator<2> generator;

        std::vector<TissueCell> cells;
        generator.GenerateForCrypt(cells, *p_mesh, false);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        double y0 = 0.3;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);

            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = cells[i].GetCellCycleModel()->GetGeneration();
            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }

            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }
    }
    
    void TestTysonNovakCellCycleModelCellsGenerator() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;

        TysonNovakCellCycleModelCellsGenerator<2> generator;

        std::vector<TissueCell> cells;
        generator.GenerateForCrypt(cells, *p_mesh, true);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);
        }
    }


    void TestWntCellCycleModelCellsGenerator() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;

        WntCellCycleModelCellsGenerator<2> generator;

        std::vector<TissueCell> cells;
        generator.GenerateForCrypt(cells, *p_mesh, false);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);

            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }
    }
    
    void TestSimpleWntCellCycleModelCellsGenerator() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;

        SimpleWntCellCycleModelCellsGenerator<2> generator;

        std::vector<TissueCell> cells;
        generator.GenerateForCrypt(cells, *p_mesh, false);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        double y0 = 0.3;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);

            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = cells[i].GetCellCycleModel()->GetGeneration();
            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }

            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }
    }
    
    
    void TestStochasticWntCellCycleModelCellsGenerator() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;

        StochasticWntCellCycleModelCellsGenerator<2> generator;

        std::vector<TissueCell> cells;
        generator.GenerateForCrypt(cells, *p_mesh, false);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        double y0 = 0.3;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetLocationIndex(), i);

            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = cells[i].GetCellCycleModel()->GetGeneration();
            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }

            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }
    }

};

#endif /*TESTCELLSGENERATOR_HPP_*/
