#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
                                   double birthTime,
                                   unsigned int generation,
                                   AbstractCellCycleModel *pCellCycleModel)
        : mpCellCycleModel(pCellCycleModel)
{
    mpSimulationTime = SimulationTime::Instance();
    // Stem cells are the only ones with generation = 0
    assert( (generation == 0) == (cellType == STEM) );
    //std::cout<< "Birth time " << birthTime << "\n" ;
    mBirthTime=birthTime;
    mGeneration=generation;
    mCellType=cellType;
    mpCellCycleModel->SetCellType(cellType);
    mCanDivide = false;
}

MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
                                   unsigned int generation,
                                   AbstractCellCycleModel *pCellCycleModel)
        : mpCellCycleModel(pCellCycleModel)
{
    mpSimulationTime = SimulationTime::Instance();
    // Stem cells are the only ones with generation = 0
    assert( (generation == 0) == (cellType == STEM) );
    
    mBirthTime=mpSimulationTime->GetDimensionalisedTime();
    mGeneration=generation;
    mCellType=cellType;
    mpCellCycleModel->SetCellType(cellType);
    mCanDivide = false;
}

void MeinekeCryptCell::CommonCopy(const MeinekeCryptCell &other_cell)
{
    // Copy 'easy' data members
    mBirthTime = other_cell.mBirthTime;
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
    mBirthTime = birthTime;
}

void MeinekeCryptCell::SetBirthTime()
{
    mBirthTime = mpSimulationTime->GetDimensionalisedTime();
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

double MeinekeCryptCell::GetAge(double simulationTime)
{
    //std::cout<< "SImulation time" << simulationTime << "\n" ;
    return simulationTime-mBirthTime;
}

double MeinekeCryptCell::GetAge()
{
    return mpSimulationTime->GetDimensionalisedTime() - mBirthTime;
}

unsigned int MeinekeCryptCell::GetGeneration()
{
    return mGeneration;
}

CryptCellType MeinekeCryptCell::GetCellType()
{
    return mCellType;
}


bool MeinekeCryptCell::ReadyToDivide()
{
    assert(mpSimulationTime!=NULL);

    //std::cout<< "Divide time" << mpSimulationTime->GetDimensionalisedTime() << "\n" ;
    mCanDivide = mpCellCycleModel->ReadyToDivide(mpSimulationTime->GetDimensionalisedTime() - mBirthTime);
    return mCanDivide;
}

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
            mBirthTime = mpSimulationTime->GetDimensionalisedTime();
            return MeinekeCryptCell(TRANSIT, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
        else
        {
            mGeneration++;
            mCellType = DIFFERENTIATED;
            mpCellCycleModel->SetCellType(mCellType);
            mBirthTime = mpSimulationTime->GetDimensionalisedTime();
            return MeinekeCryptCell(DIFFERENTIATED, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
    }
    else
    {
        mBirthTime = mpSimulationTime->GetDimensionalisedTime();
        return MeinekeCryptCell(TRANSIT, 1,
                                mpCellCycleModel->CreateCellCycleModel());
    }
    mCanDivide = false;
}
