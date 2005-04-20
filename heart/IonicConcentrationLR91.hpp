#ifndef _IONICCONCENTRATIONLR91_HPP_
#define _IONICCONCENTRATIONLR91_HPP_


class IonicConcentrationLR91
{
    protected:
        double mMagnitudeOfIonicConcentration;
    
    public:
        // Constructor
        IonicConcentrationLR91(void);
        // Destructor
        ~IonicConcentrationLR91(void);
        // Set and get methods
        void SetMagnitudeOfIonicConcentration(const double &rMagnitudeOfIonicConcentration);
        double GetMagnitudeOfIonicConcentration(void);
};


#endif //_IONICCONCENTRATION_HPP_
