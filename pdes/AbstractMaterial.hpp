#ifndef _ABSTRACTMATERIAL_HPP_
#define _ABSTRACTMATERIAL_HPP_

#if 0

class AbstractMaterial 
{
    
public:
    //    bool CheckReady();
    //    void PrintLaw();

    // need to be implemented by derived class

    virtual double GetDensity()=0;

    virtual MatrixDouble ComputeStress(MatrixDouble F);
    virtual double*      ComputeDTdE  (MatrixDouble F);

    virtual double GetZeroStrainPressure()=0;

private:
    virtual double GetdW_by_dI1   (double I1,double I2,double I3)=0;
    virtual double GetdW_by_dI2   (double I1,double I2,double I3)=0;
    virtual double GetdW_by_dI2   (double I1,double I2,double I3)=0;

    virtual double Getd2W_by_dI1  (double I1,double I2,double I3)=0;
    virtual double Getd2W_by_dI2  (double I1,double I2,double I3)=0;
    virtual double Getd2W_by_dI3  (double I1,double I2,double I3)=0;
    virtual double Getd2W_by_dI1I2(double I1,double I2,double I3)=0;
    virtual double Getd2W_by_dI2I3(double I1,double I2,double I3)=0;
    virtual double Getd2W_by_dI3I1(double I1,double I2,double I3)=0;

    virtual double GetdW_by_dE(MatrixDouble F, int index1, int index2)=0;
    
protected:

    double density;

    // can be accessed by derived class
    bool bParamSet;
    bool bDensitySet;
    bool bNonDimensionalised;
    bool bIsIncompressible;
    bool bIsTransIsoLaw;    
};


#endif //_ABSTRACTMATERIAL_HPP_
