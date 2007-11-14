#ifndef ELLIPTICPDEWITHRADIALLINEARSOURCE_HPP_
#define ELLIPTICPDEWITHRADIALLINEARSOURCE_HPP_


/**
 * The pde u_xx+u_yy - (x^2+y^2)u = 0, just for one particular test. This has the solution
 * exp(xy)
 */
class EllipticPdeWithRadialLinearSource :public AbstractLinearEllipticPde<2>
{
public:

    double ComputeConstantInUSourceTerm(ChastePoint<2> )
    {
        return 0.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(ChastePoint<2> x)
    {
        return -(x[0]*x[0] + x[1]*x[1]); 
    }
    
    c_matrix<double, 2, 2> ComputeDiffusionTerm(ChastePoint<2> )
    {
        return identity_matrix<double>(2);
    }
};

#endif /*ELLIPTICPDEWITHRADIALLINEARSOURCE_HPP_*/
