#ifndef POLEZERO3DIN1DLAW_HPP_
#define POLEZERO3DIN1DLAW_HPP_

#include "FiniteElasticityAssembler.hpp" // lazy

class PoleZero3dIn1dLaw
{
    static const double TOL = 1e-8;
    
    std::vector<double> mPStore;
    std::vector<double> mE22Store;
    std::vector<double> mE33Store;
    bool mUseStore;
    
    double mLastE22;
    double mLastE33;
    
    std::vector<double> k;
    std::vector<double> a;
    std::vector<double> b;

    double dWdE(double Eii, unsigned index)
    {
        assert(index>=0 && index<=2);
        if(Eii >= a[index])
        {
            #define COVERAGE_IGNORE
            std::stringstream ss;
            ss << "Eii > aii, where i = " << index << " and Eii = " << Eii << "\n";
            EXCEPTION(ss.str());
            #undef COVERAGE_IGNORE
        }
        
        return   k[index] 
               * Eii 
               * (2+b[index]*Eii/(a[index]-Eii))
               / pow(a[index]-Eii,b[index]);
    }

    double d2WdE2(double Eii, unsigned index)
    {
        assert(index>=0 && index<=3);
        assert(Eii < a[index]);

        return    k[index] 
                * pow(a[index]-Eii, -b[index]-2) 
                * (   2*(a[index]-Eii)*(a[index]-Eii)
                    + 4*b[index]*Eii*(a[index]-Eii)
                    + b[index]*(b[index]+1)*Eii*Eii  );
    }

    
    void CalcResidual(Tensor<1,3>& resid, const double E, double p, double E22, double E33)
    {
        resid[0] = (2*E+1)*(2*E22+1)*(2*E33+1) - 1;
        resid[1] = dWdE(E22,1) - p/(2*E22+1);
        resid[2] = dWdE(E33,2) - p/(2*E33+1);
    }
    
    double NormResidual(Tensor<1,3>& resid)
    {
        return sqrt( resid[0]*resid[0] + resid[1]*resid[1] + resid[2]*resid[2] );
    }


    void CalcJacobian(Tensor<2,3>& jac, const double E, double p, double E22, double E33)
    {
        jac[0][0] = 0;
        jac[0][1] = (2*E+1) * 2 * (2*E33+1);
        jac[0][2] = (2*E+1) * (2*E22+1) * 2;

        jac[1][0] = - 1/(2*E22+1);
        jac[1][1] = d2WdE2(E22,1) + 2*p/( (2*E22+1)*(2*E22+1) );
        jac[1][2] = 0.0;

        jac[2][0] = - 1/(2*E33+1);
        jac[2][1] = 0.0;
        jac[2][2] = d2WdE2(E33,2) + 2*p/( (2*E33+1)*(2*E33+1) );
    }

        
    double SolveSystem(const double E)
    {
        Tensor<1,3> resid;
        Tensor<2,3> jac;
        
        double p = 0.0;
        double C22 = 1/sqrt(2*E+1);
        double E22 = 0.5*(C22-1);   
        double E33 = E22;
        
        if(mUseStore && E>=-0.3)
        {
            unsigned index = (unsigned)floor( (E+0.3)*mPStore.size() );
            assert(index < mPStore.size());
            p = mPStore[index];
            E22 = mE22Store[index];
            E33 = mE33Store[index];
        }
        
        CalcResidual(resid, E, p, E22, E33);
        double norm = NormResidual(resid);
        unsigned counter = 0;

//        Tensor<2,3> num_jac;
//        double h=0.001;
//        Tensor<1,3> resid2;
//        CalcResidual(resid2, E, p+h, E22, E33);
//        for(unsigned i=0; i<3; i++)
//        {
//            num_jac[i][0] = (resid2[i]-resid[i])/h;
//        }
//        CalcResidual(resid2, E, p, E22+h, E33);
//        for(unsigned i=0; i<3; i++)
//        {
//            num_jac[i][1] = (resid2[i]-resid[i])/h;
//        }
//        CalcResidual(resid2, E, p, E22, E33+h);
//        for(unsigned i=0; i<3; i++)
//        {
//            num_jac[i][2] = (resid2[i]-resid[i])/h;
//        }
//        CalcJacobian(jac, E, p, E22, E33);
//        for(unsigned i=0; i<3; i++)
//        {
//            for(unsigned j=0; j<3; j++)
//            {
//                std::cout << i << " " << j << " " << num_jac[i][j]-jac[i][j] << ", " << num_jac[i][j] << " " << jac[i][j] << "\n";
//            }
//        }        
//assert(0);


        while ((norm > TOL) && (counter<20))
        {
            CalcJacobian(jac, E, p, E22, E33);
            Tensor<2,3> invJac = invert(jac);
            
            std::vector<double> update(3,0.0);

            for(unsigned i=0; i<3; i++)
            {
                for(unsigned j=0; j<3; j++)
                {
                   update[i] += invJac[i][j]*resid[j];
                }
            }        

            double damping = ChooseBestUpdate(E, update, p, E22, E33);

            p  -= damping*update[0];
            E22 -= damping*update[1];
            E33 -= damping*update[2];
            
            CalcResidual(resid, E, p, E22, E33);
            norm = NormResidual(resid);
            counter++;
        }

        assert(counter<20);

        mLastE22 = E22;
        mLastE33 = E33;

        return p;
    }

    double ChooseBestUpdate(double E, std::vector<double>& update, double p, double E22, double E33)
    {
        std::vector<double> try_vars = update;
        double best_damping = 0.0;

        Tensor<1,3> resid;
        CalcResidual(resid, E, p, E22, E33);
        double best_norm = NormResidual(resid);

        for(unsigned i=1; i<=10; i++)
        {
            double damping = i/10.0;

            try_vars[0] = p   - damping*update[0];
            try_vars[1] = E22 - damping*update[1];
            try_vars[2] = E33 - damping*update[2];
            
            CalcResidual(resid, E, try_vars[0], try_vars[1], try_vars[2]);
            double norm = NormResidual(resid);
            
            if(norm < best_norm)
            {
                best_norm = norm;
                best_damping = damping;
            }
        }
        
        if(best_damping == 0)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Newton step increased residual");
            #undef COVERAGE_IGNORE
        }
        
        return best_damping;
    }

public : 
    PoleZero3dIn1dLaw()
    {
        k.resize(3);
        a.resize(3);
        b.resize(3);
        
        k[0] = 2;
        k[1] = 2;
        k[2] = 2;
        a[0] = 0.476;
        a[1] = 0.619;
        a[2] = 0.943;
        b[0] = 1.5;
        b[1] = 1.5;
        b[2] = 0.442;
        
        mUseStore = false;
    }
    
    void SetUpStores()
    {
        unsigned num = 1000; 
        for(unsigned i=0; i<=num; i++)
        {
            double E =  -0.3 + i*0.3/num;
            double p = SolveSystem(E);
            
            mPStore.push_back(p);
            mE22Store.push_back(mLastE22);
            mE33Store.push_back(mLastE33);
        }
            
        mUseStore = true;       
    }
        

    double GetT(const double E)
    {
        if(E>=0)
        {
            return dWdE(E,0); 
        }
        else
        {
            double p = SolveSystem(E);
            return -p/(2*E+1);
        }
    }
    
    double GetDTdE(const double E)
    {
        double h = 0.00001;
        return (GetT(E+h)-GetT(E))/h;
    }
};



#endif /*POLEZERO3DIN1DLAW_HPP_*/
