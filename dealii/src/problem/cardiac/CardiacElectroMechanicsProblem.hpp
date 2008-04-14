/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CARDIACELECTROMECHANICSPROBLEM_HPP_
#define CARDIACELECTROMECHANICSPROBLEM_HPP_

#include "MooneyRivlinMaterialLaw.hpp"
#include "CardiacMechanicsAssembler.cpp"
#include "ImplicitCardiacMechanicsAssembler.hpp"
#include "FiniteElasticityTools.hpp"
#include "LogFile.hpp"


/**
 *  Solve a full cardiac electro-mechanics problem in 2d or 3d.
 * 
 *  See documentation for AbstractCardiacElectroMechanicsProblem
 */
template<unsigned DIM>
class CardiacElectroMechanicsProblem : public AbstractCardiacElectroMechanicsProblem<DIM>
{
private:
    unsigned mNumElementsPerDimInMechanicsMesh;

public:
    /** 
     *  Constructor
     *  @param pCellFactory cell factory for creating cells (see Monodomain tests)
     *  @endTime end time of the simulation. Start time is assumed to be 0.0
     *  @numElementsPerDimInMechanicsMesh number of elements in each direction
     *  in the mechanics mesh  
     *  @useExplicit Whether to use an explicit or implicit mechanics solver
     *  @outputDirectory. Output directory. Omit if no output is required.
     * 
     *  See documentation for AbstractCardiacElectroMechanicsProblem
     */
    CardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   unsigned numElementsPerDimInMechanicsMesh,
                                   bool useExplicitMethod,
                                   unsigned numElecStepsPerMechStep,
                                   double nhsOdeTimeStep,
                                   std::string outputDirectory = "")
        :  AbstractCardiacElectroMechanicsProblem<DIM>(pCellFactory,
                                                       endTime,
                                                       useExplicitMethod,
                                                       numElecStepsPerMechStep,
                                                       nhsOdeTimeStep,
                                                       outputDirectory)
    {
        mNumElementsPerDimInMechanicsMesh = numElementsPerDimInMechanicsMesh;
    }
    

    void ConstructMeshes()
    {        
        double width = 1.0;
        
        // create electrics mesh
        this->mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();

        unsigned num_elem = 96;
        this->mpElectricsMesh->ConstructRectangularMesh(num_elem,num_elem);
        this->mpElectricsMesh->Scale(width/num_elem,width/num_elem);

        // create mechanics mesh
        this->mpMechanicsMesh = new Triangulation<DIM>();
        Point<2> zero;
        Point<2> opposite_corner;
        opposite_corner[0] = width;
        opposite_corner[1] = width;
        
        std::vector<unsigned> repetitions;
        repetitions.push_back(mNumElementsPerDimInMechanicsMesh);
        repetitions.push_back(mNumElementsPerDimInMechanicsMesh);
        
        GridGenerator::subdivided_hyper_rectangle(*(this->mpMechanicsMesh), repetitions, zero, opposite_corner);
        
        LOG(1, "Width of meshes is " << width);
        LOG(1, "Num nodes in electrical and mechanical meshes are: " << this->mpElectricsMesh->GetNumNodes() << ", " << this->mpMechanicsMesh->n_vertices() << "\n");
    }

    
    void ConstructMechanicsAssembler(std::string mechanicsOutputDir)
    {
        Point<DIM> zero;
        FiniteElasticityTools<DIM>::SetFixedBoundary(*(this->mpMechanicsMesh), 0, 0.0); 
               
        if(this->mUseExplicitMethod)
        {
            this->mpCardiacMechAssembler = new CardiacMechanicsAssembler<DIM>(this->mpMechanicsMesh,mechanicsOutputDir);
        }
        else
        {
            ImplicitCardiacMechanicsAssembler<DIM>* p_assembler = new ImplicitCardiacMechanicsAssembler<DIM>(this->mpMechanicsMesh,mechanicsOutputDir);
            this->mpCardiacMechAssembler = p_assembler;
        }
    }
};

#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
