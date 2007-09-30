
#include "CellwiseData.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>


/** Pointer to the single instance */
template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::mpInstance = NULL;


template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CellwiseData<DIM>;
    }
    return mpInstance;
}


template<unsigned DIM>
CellwiseData<DIM>::CellwiseData()
 :  mpCrypt(NULL),
    mAllocatedMemory(false)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}

template<unsigned DIM>
CellwiseData<DIM>::~CellwiseData()
{}

template<unsigned DIM>
void CellwiseData<DIM>::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

template<unsigned DIM>
double CellwiseData<DIM>::GetValue(MeinekeCryptCell* pCell, unsigned variableNumber)
{
    assert(IsSetUp());
    assert(mpCrypt!=NULL);
    assert(mAllocatedMemory);
    assert(variableNumber < mNumberOfVariables);
    assert(pCell!=NULL);

    unsigned node_index = mpCrypt->GetNodeCorrespondingToCell(*pCell)->GetIndex();
    unsigned vector_index = node_index*mNumberOfVariables + variableNumber;
    
    return mData[vector_index];
}

template<unsigned DIM>
void CellwiseData<DIM>::SetValue(double value, Node<DIM>* pNode, unsigned variableNumber)
{
    if (mAllocatedMemory==false)
    {
        std::cout << "mAllocatedMemory = FALSE" << "\n";
    }
    if (mpInstance==NULL)
    {   
        std::cout << "mpInstance = NULL" << "\n";
    } 
    if (mpCrypt==NULL)
    {
        std::cout << "mpCrypt = NULL" << "\n";
    }
    
    assert(IsSetUp());
    assert(variableNumber < mNumberOfVariables);
    unsigned vector_index = pNode->GetIndex()*mNumberOfVariables + variableNumber;
    mData[vector_index] = value;
}

template<unsigned DIM>
void CellwiseData<DIM>::SetCrypt(Crypt<DIM>& rCrypt)
{
    if(mAllocatedMemory == false)
    {
        EXCEPTION("SetCrypt must be called after SetNumNodesAndVars()");
    }

    mpCrypt=&rCrypt;
    rCrypt.InitialiseCells();
}

template<unsigned DIM>
void CellwiseData<DIM>::SetNumNodesAndVars(unsigned numNodes, unsigned numberOfVariables)
{
    if(mpCrypt!=NULL)
    {
        EXCEPTION("SetNumNodesAndVars() must be called before setting the crypt (and after a Destroy)");
    }
    
    assert(numberOfVariables>0);
    assert(mAllocatedMemory==false);
    
    mNumberOfVariables = numberOfVariables;
    mData.clear();
    mData.resize(numNodes * mNumberOfVariables, 0.0);
    
    mAllocatedMemory = true;
}

template<unsigned DIM>
bool CellwiseData<DIM>::IsSetUp()
{
    return ((mAllocatedMemory) && (mpInstance!=NULL) && (mpCrypt!=NULL));   
}

