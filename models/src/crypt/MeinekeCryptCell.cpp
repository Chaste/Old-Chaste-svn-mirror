#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


//MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
//                                   double birthTime,
//                                   unsigned int generation,
//                                   AbstractCellCycleModel *pCellCycleModel)
//        : mpCellCycleModel(pCellCycleModel)
//{
//    mpSimulationTime = SimulationTime::Instance();
//    // Stem cells are the only ones with generation = 0
//    assert( (generation == 0) == (cellType == STEM) );
//    //std::cout<< "Birth time " << birthTime << "\n" ;
//    mBirthTime=birthTime;
//    mGeneration=generation;
//    mCellType=cellType;
//    mpCellCycleModel->SetCellType(cellType);
//    mCanDivide = false;
//}

MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
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
    assert( (generation == 0) == (cellType == STEM) );
    mGeneration=generation;
    mCellType=cellType;
    mpCellCycleModel->SetCellType(cellType);
    mCanDivide = false;
}

void MeinekeCryptCell::CommonCopy(const MeinekeCryptCell &other_cell)
{
    // Copy 'easy' data members
    mGeneration = other_cell.mGeneration;
    mCellType = other_cell.mCellType;
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

unsigned int MeinekeCryptCell::GetGeneration()
{
    return mGeneration;
}

CryptCellType MeinekeCryptCell::GetCellType()
{
    return mCellType;
}


bool MeinekeCryptCell::ReadyToDivide(std::vector<double> cellCycleInfluences)
{	
	mCanDivide = mpCellCycleModel->ReadyToDivide(cellCycleInfluences);
    return mCanDivide;
}

//bool MeinekeCryptCell::ReadyToDivide(double wnt_stimulus)
//{
//    mCanDivide = mpCellCycleModel->ReadyToDivide(wnt_stimulus);
//    return mCanDivide;
//}

MeinekeCryptCell MeinekeCryptCell::Divide()
{
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
            return MeinekeCryptCell(TRANSIT, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
        else
        {
            mGeneration++;
            mCellType = DIFFERENTIATED;
            mpCellCycleModel->SetCellType(mCellType);
            mpCellCycleModel->ResetModel();// Cell goes back to age zero
            return MeinekeCryptCell(DIFFERENTIATED, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
    }
    else
    {
        mpCellCycleModel->ResetModel();// Cell goes back to age zero
        return MeinekeCryptCell(TRANSIT, 1,
                                mpCellCycleModel->CreateCellCycleModel());
    }
    mCanDivide = false;
}
