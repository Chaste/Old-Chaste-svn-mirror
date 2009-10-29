
#include "ExplicitCardiacMechanicsAssembler.hpp"

template<unsigned DIM>
ExplicitCardiacMechanicsAssembler<DIM>::ExplicitCardiacMechanicsAssembler(ContractionModel contractionModel,
                                                                          QuadraticMesh<DIM>* pQuadMesh,
                                                                          std::string outputDirectory,
                                                                          std::vector<unsigned>& rFixedNodes,
                                                                          AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
    : AbstractCardiacMechanicsAssembler<DIM>(pQuadMesh,
                                             outputDirectory,
                                             rFixedNodes,
                                             pMaterialLaw)
{
    switch(contractionModel)
    {
        case TEST1:
        {
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                this->mContractionModelSystems.push_back(new NonPhysiologicalContractionModel(1));
            }
            break;
        }
        case KERCHOFFS2003: //stretch dependent, will this work with explicit??
        {
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                Kerchoffs2003ContractionModel* p_model = new Kerchoffs2003ContractionModel();
                this->mContractionModelSystems.push_back(p_model);
            }
            break;
        }
        default:
        {
            EXCEPTION("Unknown or stretch-rate-dependent contraction model");
        } 
    }
       
    assert(!(this->mContractionModelSystems[0]->IsStretchRateDependent()));
}

template<unsigned DIM>
ExplicitCardiacMechanicsAssembler<DIM>::~ExplicitCardiacMechanicsAssembler()
{        
    for(unsigned i=0; i<this->mContractionModelSystems.size(); i++)
    {
        //// memory leak as this is commented out. But get glibc failure with it in... (EMTODO2)
        //delete this->mContractionModelSystems[i];
    }
}        
    
template<unsigned DIM>
void ExplicitCardiacMechanicsAssembler<DIM>::GetActiveTensionAndTensionDerivs(double currentFibreStretch, 
                                                                              unsigned currentQuadPointGlobalIndex,
                                                                              bool assembleJacobian,
                                                                              double& rActiveTension,
                                                                              double& rDerivActiveTensionWrtLambda,
                                                                              double& rDerivActiveTensionWrtDLambdaDt)
{
    // the active tensions have already been computed for each contraction model, so can 
    // return it straightaway..
    rActiveTension = this->mContractionModelSystems[currentQuadPointGlobalIndex]->GetActiveTension();

    // these are unset
    rDerivActiveTensionWrtLambda = 0.0;
    rDerivActiveTensionWrtDLambdaDt = 0.0;

    // store the value of given for this quad point, so that it can be used when computing 
    // the active tension at the next timestep
    this->mStretches[currentQuadPointGlobalIndex] = currentFibreStretch;
}

template<unsigned DIM>
void ExplicitCardiacMechanicsAssembler<DIM>::Solve(double time, double nextTime, double odeTimestep)
{
    assert(time < nextTime);
    this->mCurrentTime = time;
    this->mNextTime = nextTime;
    this->mOdeTimestep = odeTimestep;        

    // assemble the residual again so that mStretches is set (in GetActiveTensionAndTensionDerivs)
    // using the current deformation.
    this->AssembleSystem(true,false);
            
    // integrate contraction models
    for(unsigned i=0; i<this->mContractionModelSystems.size(); i++)
    {
        this->mContractionModelSystems[i]->SetStretchAndStretchRate(this->mStretches[i], 0.0);
        this->mContractionModelSystems[i]->RunDoNotUpdate(time, nextTime, odeTimestep);
        this->mContractionModelSystems[i]->UpdateStateVariables();
    }   
    
    // solve
    NonlinearElasticityAssembler<DIM>::Solve();
}



template class ExplicitCardiacMechanicsAssembler<2>;
template class ExplicitCardiacMechanicsAssembler<3>;
