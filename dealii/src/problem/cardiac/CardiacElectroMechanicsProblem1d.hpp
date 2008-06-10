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


#ifndef CARDIACELECTROMECHANICSPROBLEM1D_HPP_
#define CARDIACELECTROMECHANICSPROBLEM1D_HPP_

#include "AbstractCardiacElectroMechanicsProblem.hpp"
#include "Explicit1dCardiacMechanicsAssembler.hpp"
#include "Implicit1dCardiacMechanicsAssembler.hpp"



/**
 *  A 1d CardiacElectroMechanics assembler
 *
 *  Note 1d incompressible mechanics doesn't any sense, we can't just
 *  use CardiacMechanicsAssembler<1>. Instead a special 1d cardiac mechanics
 *  assembler, which uses a particular material law that takes uni-axial
 *  deformation in 3d and returns the corresponding 1d stress, is used. An
 *  implicit or explicit version can be used.
 *
 *  See also AbstractCardiacElectroMechanicsProblem
 */
class CardiacElectroMechanicsProblem1d : public AbstractCardiacElectroMechanicsProblem<1>
{
private:
    out_stream mpFibreLengthFile;

    /** Overloaded PostSolve() writing the length of the fibre at each time to a
     *  file.
     */
    void PostSolve(double currentTime)
    {
        if(!(this->mWriteOutput))
        {
            return;
        }

        std::vector<Vector<double> >& r_deformed_solution
         = dynamic_cast<AbstractElasticityAssembler<1>*>
           (this->mpCardiacMechAssembler)->rGetDeformedPosition();

        assert(r_deformed_solution.size()==1);

        double length = -1;
        for(unsigned i=0; i<r_deformed_solution[0].size(); i++)
        {
            if(r_deformed_solution[0](i)>length)
            {
                length = r_deformed_solution[0](i);
            }
        }

        // verify we found something
        assert(length>0);

        mpFibreLengthFile->precision(8);
        (*mpFibreLengthFile) << currentTime << " " << length << "\n";
    }


public:
    CardiacElectroMechanicsProblem1d(AbstractCardiacCellFactory<1>* pCellFactory,
                                     double endTime,
                                     unsigned numElecStepsPerMechStep,
                                     std::string outputDirectory = "")
        :  AbstractCardiacElectroMechanicsProblem<1>(pCellFactory,
                                                     endTime,
                                                     numElecStepsPerMechStep,
                                                     0.01,
                                                     outputDirectory)
    {
        if(this->mWriteOutput)
        {
            OutputFileHandler output_file_handler(this->mOutputDirectory, false);
            mpFibreLengthFile = output_file_handler.OpenOutputFile("length.txt");
        }
    }

    ~CardiacElectroMechanicsProblem1d()
    {
        if(this->mWriteOutput)
        {
            mpFibreLengthFile->close();
        }
    }

    void ConstructMeshes()
    {
        // create electrics mesh
        mpElectricsMesh = new ConformingTetrahedralMesh<1,1>();
        unsigned num_elem = 128;

        mpElectricsMesh->ConstructLinearMesh(num_elem);
        mpElectricsMesh->Scale(1.0/num_elem);

        // create mechanics mesh
        mpMechanicsMesh = new Triangulation<1>();
        GridGenerator::hyper_cube(*mpMechanicsMesh, 0.0, 1.0);
        mpMechanicsMesh->refine_global(7);

        std::cout << "Number of nodes = " << mpElectricsMesh->GetNumNodes() << ", " << mpMechanicsMesh->n_vertices() << "\n";

        assert(mpMechanicsMesh->n_vertices()==mpElectricsMesh->GetNumNodes());
    }


    void ConstructMechanicsAssembler(std::string mechanicsOutputDir)
    {
        mpCardiacMechAssembler = new Implicit1dCardiacMechanicsAssembler(mpMechanicsMesh, mechanicsOutputDir);
    }
};


#endif /*CARDIACELECTROMECHANICSPROBLEM1D_HPP_*/
