#include "MeinekeCryptCell.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
								   CryptCellMutationState mutationState,
                                   unsigned int generation,
                                   AbstractCellCycleModel *pCellCycleModel)
        : mpCellCycleModel(pCellCycleModel)
{
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsSimulationTimeSetUp()==false)
    {
    	EXCEPTION("MeinekeCryptCell is setting up a cell cycle model but SimulationTime has not been set up");	
    }
    // Stem cells are the only ones with generation = 0
    //assert( (generation == 0) == (cellType == STEM) ); Not for Wnt cells
    mGeneration=generation;
    mCellType=cellType;
    mMutationState = mutationState;
    mpCellCycleModel->SetCellType(cellType);
    mCanDivide = false;
}

void MeinekeCryptCell::CommonCopy(const MeinekeCryptCell &other_cell)
{
	// Copy 'easy' data members
	mGeneration = other_cell.mGeneration;
    mCellType = other_cell.mCellType;
    mMutationState = other_cell.mMutationState;
    mCanDivide = other_cell.mCanDivide;
    mpSimulationTime = other_cell.mpSimulationTime;
    // Copy cell cycle model
    // First create a new object
    mpCellCycleModel = other_cell.mpCellCycleModel->CreateCellCycleModel();
    // Then copy its state
    *mpCellCycleModel = *(other_cell.mpCellCycleModel);
    
}

MeinekeCryptCell::MeinekeCryptCell(const MeinekeCryptCell &other_cell)
{
    CommonCopy(other_cell);
}

void MeinekeCryptCell::operator=(const MeinekeCryptCell &other_cell)
{
    delete mpCellCycleModel;
    CommonCopy(other_cell);
}

MeinekeCryptCell::~MeinekeCryptCell()
{
    delete mpCellCycleModel;
}


void MeinekeCryptCell::SetBirthTime(double birthTime)
{
    mpCellCycleModel->SetBirthTime(birthTime);
}

void MeinekeCryptCell::SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel)
{
    delete mpCellCycleModel;
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCellType(mCellType);
}

AbstractCellCycleModel *MeinekeCryptCell::GetCellCycleModel()
{
    return mpCellCycleModel;
}

void MeinekeCryptCell::SetNodeIndex(unsigned index)
{
    mNodeIndex = index;
}

unsigned MeinekeCryptCell::GetNodeIndex()
{
    return mNodeIndex;
}

double MeinekeCryptCell::GetAge()
{
	return mpCellCycleModel->GetAge();
}

double MeinekeCryptCell::GetBirthTime()
{
	return mpCellCycleModel->GetBirthTime();
}

unsigned int MeinekeCryptCell::GetGeneration()
{
    return mGeneration;
}

CryptCellType MeinekeCryptCell::GetCellType()
{
    return mCellType;
}

CryptCellMutationState MeinekeCryptCell::GetMutationState()
{
    return mMutationState;
}

void MeinekeCryptCell::SetCellType(CryptCellType cellType)
{
	mCellType = cellType;
	mpCellCycleModel->SetCellType(mCellType);
}

void MeinekeCryptCell::SetMutationState(CryptCellMutationState mutationState)
{
	mMutationState = mutationState;
}

/**
 * The MeinekeCryptCell ready to divide method
 * 
 * @param cellCycleInfluences a std::vector of doubles, with any relevant cell 
 * cycle influences in it. This function pushes back the mutation state of this 
 * cell onto the vector before sending it to the cell cycle models.
 */
bool MeinekeCryptCell::ReadyToDivide(std::vector<double> cellCycleInfluences)
{	
	double mutation_state = -1;
	if(mMutationState==HEALTHY)
	{
		//std::cout << "HEALTHY" << std::endl;
		mutation_state=0;	
	}
	if(mMutationState==APC_ONE_HIT)	
	{
		//std::cout << "APC +/-" << std::endl;
		mutation_state=1;	
	}
	if(mMutationState==BETA_CATENIN_ONE_HIT)	
	{
		//std::cout << "Beta-cat +/-" << std::endl;
		mutation_state=2;	
	}
	if(mMutationState==APC_TWO_HIT)	
	{
		//std::cout << "APC -/-" << std::endl;
		mutation_state=3;	
	}
	if(fabs(mutation_state+1)<1e-6)
	{
#define COVERAGE_IGNORE
		EXCEPTION("This cell has an invalid mutation state");
#undef COVERAGE_IGNORE
	}
	cellCycleInfluences.push_back(mutation_state);
	mCanDivide = mpCellCycleModel->ReadyToDivide(cellCycleInfluences);
    return mCanDivide;
}

MeinekeCryptCell MeinekeCryptCell::Divide()
{	//Copy this cell and give new one relevant attributes...
	assert(mpSimulationTime!=NULL);
    assert(mCanDivide);
    CancerParameters *p_params = CancerParameters::Instance();
    //std::cout<< "Divide time" << mpSimulationTime->GetDimensionalisedTime() << "\n" ;
    if (mCellType != STEM)
    {
        if (mGeneration < p_params->GetMaxTransitGenerations())
        {
            mGeneration++;
            mpCellCycleModel->ResetModel();// Cell goes back to age zero
            return MeinekeCryptCell(TRANSIT, mMutationState, mGeneration, 
                                    mpCellCycleModel->CreateCellCycleModel());
        }
        else
        {
            mGeneration++;
            mCellType = DIFFERENTIATED;
            mpCellCycleModel->SetCellType(mCellType);
            mpCellCycleModel->ResetModel();// Cell goes back to age zero
            return MeinekeCryptCell(DIFFERENTIATED, mMutationState, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
    }
    else
    {
        mpCellCycleModel->ResetModel();// Cell goes back to age zero
        return MeinekeCryptCell(TRANSIT, mMutationState, 1, 
                                mpCellCycleModel->CreateCellCycleModel());
    }
    mCanDivide = false;
}
