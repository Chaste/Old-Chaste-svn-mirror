#ifndef CARDIACELECTROMECHANICSPROBLEM1D_HPP_
#define CARDIACELECTROMECHANICSPROBLEM1D_HPP_

#include "AbstractCardiacElectroMechanicsProblem.hpp"
#include "Explicit1dCardiacMechanicsAssembler.hpp"
#include "Implicit1dCardiacMechanicsAssembler.hpp"


class CardiacElectroMechanicsProblem1d : public AbstractCardiacElectroMechanicsProblem<1>
{
public:
    CardiacElectroMechanicsProblem1d(AbstractCardiacCellFactory<1>* pCellFactory,
                                     double endTime,
                                     double timeStep,
                                     bool useExplicitMethod,
                                     std::string outputDirectory = "")
        :  AbstractCardiacElectroMechanicsProblem<1>(pCellFactory,
                                                     endTime,
                                                     timeStep,
                                                     useExplicitMethod,
                                                     outputDirectory)
    {
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
        
        std::cout << "numnodes = " << mpElectricsMesh->GetNumNodes() << ", " << mpMechanicsMesh->n_vertices() << "\n";
        
        assert(mpMechanicsMesh->n_vertices()==mpElectricsMesh->GetNumNodes());
    }
    
    
    void ConstructMechanicsAssembler()
    {
        if(mUseExplicitMethod)
        {
            mpCardiacMechAssembler = new Explicit1dCardiacMechanicsAssembler(mpMechanicsMesh);
        }
        else
        {
            mpCardiacMechAssembler = new Implicit1dCardiacMechanicsAssembler(mpMechanicsMesh);
        }
    }
};


#endif /*CARDIACELECTROMECHANICSPROBLEM1D_HPP_*/
