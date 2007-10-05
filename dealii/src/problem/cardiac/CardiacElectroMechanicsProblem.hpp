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
    unsigned mNumElementsPerDimInElectricsMesh;
    unsigned mNumElementsPerDimInMechanicsMesh;

public:
    /** 
     *  Constructor
     *  @param pCellFactory cell factory for creating cells (see Monodomain tests)
     *  @endTime end time of the simulation. Start time is assumed to be 0.0
     *  @timeStep time step for the electrics (and currently the mechanics too)
     *  @useExplicit Whether to use an explicit or implicit mechanics solver
     *  @outputDirectory. Output directory. Omit if no output is required.
     * 
     *  See documentation for AbstractCardiacElectroMechanicsProblem
     */
    CardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   double timeStep,
                                   bool useExplicitMethod,
                                   unsigned numElementsPerDimInElectricsMesh,
                                   unsigned numElementsPerDimInMechanicsMesh,
                                   std::string outputDirectory = "")
        :  AbstractCardiacElectroMechanicsProblem<DIM>(pCellFactory,
                                                       endTime,
                                                       timeStep,
                                                       useExplicitMethod,
                                                       outputDirectory)
    {
        assert(numElementsPerDimInElectricsMesh > 8);
        assert(numElementsPerDimInMechanicsMesh > 4);

        mNumElementsPerDimInElectricsMesh = numElementsPerDimInElectricsMesh;
        mNumElementsPerDimInMechanicsMesh = numElementsPerDimInMechanicsMesh;
    }
    

    void ConstructMeshes()
    {        
        double width = 0.1;
        
        // create electrics mesh
        this->mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();

        unsigned num_elem = 16; //mNumElementsPerDimInElectricsMesh;
        this->mpElectricsMesh->ConstructRectangularMesh(num_elem,num_elem);
        this->mpElectricsMesh->Scale(width/num_elem,width/num_elem);

        // create mechanics mesh
        this->mpMechanicsMesh = new Triangulation<DIM>();
        GridGenerator::hyper_cube(*(this->mpMechanicsMesh), 0.0, width);
        this->mpMechanicsMesh->refine_global(4);
        
        LOG(1, "Width of meshes is " << width);
        LOG(1, "Num nodes in electrical and mechanical meshes are: " << this->mpElectricsMesh->GetNumNodes() << ", " << this->mpMechanicsMesh->n_vertices() << "\n");
    }

    
    void ConstructMechanicsAssembler(std::string mechanicsOutputDir)
    {
        Point<DIM> zero;
        FiniteElasticityTools<DIM>::FixFacesContainingPoint(*(this->mpMechanicsMesh), zero);
        
        MooneyRivlinMaterialLaw<DIM>* p_law = new MooneyRivlinMaterialLaw<DIM>(20.0);
        if(this->mUseExplicitMethod)
        {
            this->mpCardiacMechAssembler = new CardiacMechanicsAssembler<DIM>(this->mpMechanicsMesh,mechanicsOutputDir,p_law);
        }
        else
        {
            this->mpCardiacMechAssembler = new ImplicitCardiacMechanicsAssembler<DIM>(this->mpMechanicsMesh,mechanicsOutputDir,p_law);
        }
    }
};

#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
