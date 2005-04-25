#include "CardiacCellModel.hpp"

/**
 * Constructor
 */
CardiacCellModel::CardiacCellModel()
{	
	mTransmembranePotential = 0.0;
	mStimulusCurrent = 0.0;
	mTotalIonicCurrent = 0.0;
}

/**
 * Overloaded Constructor
 * 
 * @param transmembranePotential Initial transmembrane potential
 */
CardiacCellModel::CardiacCellModel(double transmembranePotential)
{	
	mTransmembranePotential = transmembranePotential;
	mStimulusCurrent = 0.0;
	mTotalIonicCurrent = 0.0;
}

/**
 * Destructor
 */
CardiacCellModel::~CardiacCellModel()
{	
	// Do nothing
}

/**
 * Get transmembrane potential
 * 
 * @return double Current transmembrane potential
 */
double CardiacCellModel::GetTransmembranePotential(void)
{
	return mTransmembranePotential;
}

/**
 * Get stimulus current
 * 
 * @return double Stimulus current
 */
double CardiacCellModel::GetStimulusCurrent(void)
{
	return mStimulusCurrent;
}

/**
 * Get total ionic current
 * 
 * @return double Total ionic current crossing the membrane
 */
double CardiacCellModel::GetTotalIonicCurrent(void)
{
	return mTotalIonicCurrent;
}

/**
 * Sets stimulus current
 * 
 * @param stimulusCurrent New stimulus current
 */
void CardiacCellModel::SetStimulusCurrent(double stimulusCurrent)
{
	mStimulusCurrent = stimulusCurrent;
}
