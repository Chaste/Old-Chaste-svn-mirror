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
                                   unsigned numElementsPerDimInMechanicsMesh,
                                   std::string outputDirectory = "")
        :  AbstractCardiacElectroMechanicsProblem<DIM>(pCellFactory,
                                                       endTime,
                                                       timeStep,
                                                       useExplicitMethod,
                                                       outputDirectory)
    {
        //mNumElementsPerDimInElectricsMesh = numElementsPerDimInElectricsMesh;
        mNumElementsPerDimInMechanicsMesh = numElementsPerDimInMechanicsMesh;
    }
    

    void ConstructMeshes()
    {        
        double width = 1;
        
        // create electrics mesh
        this->mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();

        unsigned num_elem = 96; //mNumElementsPerDimInElectricsMesh;
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
            p_assembler->SetScaling(10);
            this->mpCardiacMechAssembler = p_assembler;
        }
    }
};

#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
