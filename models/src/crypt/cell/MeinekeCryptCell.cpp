#include "MeinekeCryptCell.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
                                   CryptCellMutationState mutationState,
                                   unsigned generation,
                                   AbstractCellCycleModel *pCellCycleModel)
        : mpCellCycleModel(pCellCycleModel)
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
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
    mUndergoingApoptosis = false;
    mIsDead = false;
    mDeathTime = DBL_MAX; // this has to be initialised for archiving...
}

void MeinekeCryptCell::CommonCopy(const MeinekeCryptCell &other_cell)
{
    // Copy 'easy' data members
    mGeneration = other_cell.mGeneration;
    mCellType = other_cell.mCellType;
    mMutationState = other_cell.mMutationState;
    mCanDivide = other_cell.mCanDivide;
    mUndergoingApoptosis = other_cell.mUndergoingApoptosis;
    mIsDead = other_cell.mIsDead;

    mDeathTime = other_cell.mDeathTime;
    mNodeIndex = other_cell.mNodeIndex;
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

MeinekeCryptCell& MeinekeCryptCell::operator=(const MeinekeCryptCell &other_cell)
{
    // In case this is self-assignment, don't delete the cell cycle model...
    AbstractCellCycleModel* temp = mpCellCycleModel;
    CommonCopy(other_cell);
    // ...until after we've copied it.
    delete temp;
    return *this;
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
    if (mpCellCycleModel != pCellCycleModel)
    {
        delete mpCellCycleModel;
    }
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCellType(mCellType);
}

AbstractCellCycleModel *MeinekeCryptCell::GetCellCycleModel() const
{
    return mpCellCycleModel;
}

/**
 * Set the node at which this cell is positioned.
 * 
 * @param index Index of the node
 */

void MeinekeCryptCell::SetNodeIndex(unsigned index)
{
    mNodeIndex = index;
}

unsigned MeinekeCryptCell::GetNodeIndex() const
{
    return mNodeIndex;
}

double MeinekeCryptCell::GetAge() const
{
    return mpCellCycleModel->GetAge();
}

double MeinekeCryptCell::GetBirthTime() const
{
    return mpCellCycleModel->GetBirthTime();
}

unsigned MeinekeCryptCell::GetGeneration() const
{
    return mGeneration;
}

CryptCellType MeinekeCryptCell::GetCellType() const
{
    return mCellType;
}

CryptCellMutationState MeinekeCryptCell::GetMutationState() const
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
    assert(!IsDead());
    if (mUndergoingApoptosis)
    {
        return false;
    }    
    
    double mutation_state = -1;
    if (mMutationState==HEALTHY || mMutationState==LABELLED)
    {
        //std::cout << "HEALTHY" << std::endl;
        mutation_state=0;
    }
    if (mMutationState==APC_ONE_HIT)
    {
        //std::cout << "APC +/-" << std::endl;
        mutation_state=1;
    }
    if (mMutationState==BETA_CATENIN_ONE_HIT)
    {
        //std::cout << "Beta-cat +/-" << std::endl;
        mutation_state=2;
    }
    if (mMutationState==APC_TWO_HIT)
    {
        //std::cout << "APC -/-" << std::endl;
        mutation_state=3;
    }
    if (fabs(mutation_state+1)<1e-6)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("This cell has an invalid mutation state");
        #undef COVERAGE_IGNORE
    }
    cellCycleInfluences.push_back(mutation_state);
    
    mCanDivide = mpCellCycleModel->ReadyToDivide(cellCycleInfluences);

    return mCanDivide;
}


void MeinekeCryptCell::StartApoptosis()
{
    assert(!IsDead());
    
    if (mUndergoingApoptosis)
    {
        EXCEPTION("StartApoptosis() called when already undergoing apoptosis");
    }
    mUndergoingApoptosis = true;
    
    CancerParameters *p_params = CancerParameters::Instance();
    
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    
    mDeathTime = p_simulation_time->GetDimensionalisedTime() + p_params->GetApoptosisTime();
}


bool MeinekeCryptCell::HasApoptosisBegun() const
{
    return mUndergoingApoptosis;
}

double MeinekeCryptCell::TimeUntilDeath() const
{
    if (!mUndergoingApoptosis)
    {
        EXCEPTION("Shouldn't be checking time until apoptosis as it isn't undergoing apoptosis");
    }
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    return mDeathTime - p_simulation_time->GetDimensionalisedTime();
}

bool MeinekeCryptCell::IsDead() const
{
    SimulationTime *p_simulation_time = SimulationTime::Instance();

    return ( mIsDead || ( (mUndergoingApoptosis) && (p_simulation_time->GetDimensionalisedTime() >= mDeathTime)) );
}

void MeinekeCryptCell::Kill()
{
    mIsDead = true;
}


MeinekeCryptCell MeinekeCryptCell::Divide()
{
    assert(!IsDead());
    
    //Copy this cell and give new one relevant attributes...
    assert(mCanDivide);
    mCanDivide = false;
    
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
    
}


void MeinekeCryptCell::UpdateCellType()
{
    mCellType = mpCellCycleModel->UpdateCellType();
}

