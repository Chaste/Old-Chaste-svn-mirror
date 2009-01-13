#ifndef CARDIACELECTROMECHPROBREGULARGEOM_HPP_
#define CARDIACELECTROMECHPROBREGULARGEOM_HPP_

#include "CardiacElectroMechanicsProblem.hpp"

template<unsigned DIM>
class CardiacElectroMechProbRegularGeom : public CardiacElectroMechanicsProblem<DIM>
{
public:
    CardiacElectroMechProbRegularGeom(double width, 
                                      unsigned numMechanicsElementsEachDir,
                                      unsigned numElectricsElementsEachDir,
                                      AbstractCardiacCellFactory<DIM>* pCellFactory,
                                      double endTime,
                                      unsigned numElecTimeStepsPerMechTimestep,
                                      double nhsOdeTimeStep,
                                      std::string outputDirectory = "")
        : CardiacElectroMechanicsProblem<DIM>(NULL, NULL, pCellFactory, endTime, 
                                              numElecTimeStepsPerMechTimestep,
                                              nhsOdeTimeStep, outputDirectory)
    {
        assert(DIM==2); // the below assumes DIM==2 

        assert(width > 0.0);
        assert(numMechanicsElementsEachDir > 0);
        assert(numElectricsElementsEachDir > 0);

        // create electrics mesh
        this->mpElectricsMesh = new TetrahedralMesh<DIM,DIM>();

        this->mpElectricsMesh->ConstructRectangularMesh(numElectricsElementsEachDir,numElectricsElementsEachDir);
        this->mpElectricsMesh->Scale(width/numElectricsElementsEachDir,width/numElectricsElementsEachDir);

        // create mechanics mesh
        this->mpMechanicsMesh = new QuadraticMesh<DIM>(width,width,numMechanicsElementsEachDir,numMechanicsElementsEachDir);
        LOG(2, "Width of meshes is " << width);
        LOG(2, "Num nodes in electrical and mechanical meshes are: " << this->mpElectricsMesh->GetNumNodes() << ", " << this->mpMechanicsMesh->GetNumNodes() << "\n");
    }
};

#endif /*CARDIACELECTROMECHPROBREGULARGEOM_HPP_*/
