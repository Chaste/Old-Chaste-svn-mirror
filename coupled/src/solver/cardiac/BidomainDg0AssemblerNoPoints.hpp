#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "BidomainPdeNoPoints.hpp"
//#include "AbstractBasisFunction.hpp"
//#include "GaussianQuadratureRule.hpp"
//#include "AbstractLinearDynamicProblemAssembler.hpp"

#include "BoundaryConditionsContainer.hpp"

#include "SimpleLinearSolver.cpp"
//#include "LinearBasisFunction.cpp"
#include "ReplicatableVector.hpp"

#define PROBLEM_DIM 2

#include <cmath>
#include "Point.hpp"
#include "Exception.hpp"
#include "UblasCustomFunctions.hpp"

/**
 * This class encapsulates tables of gaussian quadrature points and the
 * associated weights.
 *
 * Data is available for 1d, 2d and 3d quadrature over (canonical) triangles,
 * with between 1 and 3 (inclusive) gauss points in each dimension.
 */

template<int ELEM_DIM>
class GaussianQuadratureRule
{
    int mNumQuadPoints;
    std::vector<double>            mWeights;
    std::vector<c_vector<double, ELEM_DIM> >  mPoints;
    
    c_vector<double, ELEM_DIM> MakePoint(double x=0, double y=0, double z=0)
    {
        c_vector<double, ELEM_DIM> point;
        
        switch (ELEM_DIM)
        {
            case 3: point[2] = z;
            case 2: point[1] = y;
            case 1: point[0] = x;
        }
        return point;
    }
    
    
public:

    /**
     * The constructor builds the appropriate table for the dimension (given
     * by the template argument) and number of points in each dimension (given
     * as a constructor argument).
     * 
     * An exception is thrown if data is not available for the requested
     * parameters.
     */
    GaussianQuadratureRule(int numPointsInEachDimension)
    {
        mNumQuadPoints = (int) pow((double) numPointsInEachDimension,(ELEM_DIM));
        
        mWeights.reserve(mNumQuadPoints);
        mPoints.reserve(mNumQuadPoints);
        
        switch (ELEM_DIM)
        {
            case 0 :
            {
                // mNumQuadPoints should have been set to be one as
                // it is numPointsInEachDim^{0}
                mWeights.push_back(1);
                //mPoints.push_back(c_vector<double, ELEM_DIM>);
            }
            break;
            case 1 :
            {
                switch (numPointsInEachDimension)
                {
                    case 1: // 1d, 1 point
                        mWeights.push_back(1);
                        mPoints.push_back(MakePoint(0.5)); //check
                        break;
                        
                    case 2: // 1d, 2 points
                        mWeights.push_back(0.5);
                        mWeights.push_back(0.5);
                        
                        mPoints.push_back(MakePoint(0.21132486540519));
                        mPoints.push_back(MakePoint(0.78867513459481));
                        break;
                        
                    case 3: // 1d, 3 points
                        mWeights.push_back(5.0/18.0);
                        mWeights.push_back(4.0/9.0);
                        mWeights.push_back(5.0/18.0);
                        
                        mPoints.push_back(MakePoint(0.1127016654));
                        mPoints.push_back(MakePoint(0.5));
                        mPoints.push_back(MakePoint(0.8872983346));
                        break;
                        
                    default:
                        EXCEPTION("Number of gauss points per dimension not supported.");
                }
            }
            break;
            case 2 :
            {
                switch (numPointsInEachDimension)
                {
                    case 1: // 2d, 1 point per dimension
                        mWeights.push_back(0.5);
                        mPoints.push_back(MakePoint(0.25,0.5));
                        break;
                        
                    case 2: // 2d, 2 points per dimension
                        mWeights.push_back(0.19716878364870);
                        mWeights.push_back(0.19716878364870);
                        mWeights.push_back(0.05283121635130);
                        mWeights.push_back(0.05283121635130);
                        
                        mPoints.push_back(MakePoint(0.16666666666667,0.21132486540519));
                        mPoints.push_back(MakePoint(0.62200846792815,0.21132486540519));
                        mPoints.push_back(MakePoint(0.04465819873852,0.78867513459481));
                        mPoints.push_back(MakePoint(0.16666666666667,0.78867513459481));
                        break;
                        
                    case 3: // 2d, 3 points per dimension
                        mWeights.push_back(0.06846437766975);
                        mWeights.push_back(0.10954300427160);
                        mWeights.push_back(0.06846437766975);
                        mWeights.push_back(0.06172839506173);
                        mWeights.push_back(0.09876543209877);
                        mWeights.push_back(0.06172839506173);
                        mWeights.push_back(0.00869611615741);
                        mWeights.push_back(0.01391378585185);
                        mWeights.push_back(0.00869611615741);
                        
                        mPoints.push_back(MakePoint(0.10000000001607,0.11270166540000));
                        mPoints.push_back(MakePoint(0.44364916730000,0.11270166540000));
                        mPoints.push_back(MakePoint(0.78729833458393,0.11270166540000));
                        mPoints.push_back(MakePoint(0.05635083270000,0.50000000000000));
                        mPoints.push_back(MakePoint(0.25000000000000,0.50000000000000));
                        mPoints.push_back(MakePoint(0.44364916730000,0.50000000000000));
                        mPoints.push_back(MakePoint(0.01270166538393,0.88729833460000));
                        mPoints.push_back(MakePoint(0.05635083270000,0.88729833460000));
                        mPoints.push_back(MakePoint(0.10000000001607,0.88729833460000));
                        break;
                        
                    default:
                        EXCEPTION("Number of gauss points per dimension not supported.");
                }
                
            }
            break;
            case 3 :
            {
                switch (numPointsInEachDimension)
                {
                    case 1: //3d, 1 point per dimension
                        mWeights.push_back(0.12500000000000);
                        mPoints.push_back(MakePoint(0.25000000000000,0.50000000000000,0.12500000000000));
                        break;
                        
                    case 2: //3d, 2 points per dimension
                        mWeights.push_back(0.06132032652029);
                        mWeights.push_back(0.01643073197073);
                        mWeights.push_back(0.00440260136261);
                        mWeights.push_back(0.00117967347971);
                        mWeights.push_back(0.06132032652029);
                        mWeights.push_back(0.01643073197073);
                        mWeights.push_back(0.00440260136261);
                        mWeights.push_back(0.00117967347971);
                        
                        mPoints.push_back(MakePoint(0.16666666666667,   0.21132486540519,   0.13144585576580));
                        mPoints.push_back(MakePoint(0.62200846792815,   0.21132486540519,   0.03522081090086));
                        mPoints.push_back(MakePoint(0.04465819873852,   0.78867513459481,   0.03522081090086));
                        mPoints.push_back(MakePoint(0.16666666666667,   0.78867513459481,   0.00943738783766));
                        mPoints.push_back(MakePoint(0.16666666666667,   0.21132486540519,   0.49056261216234));
                        mPoints.push_back(MakePoint(0.62200846792815,   0.21132486540519,   0.13144585576580));
                        mPoints.push_back(MakePoint(0.04465819873852,   0.78867513459481,   0.13144585576580));
                        mPoints.push_back(MakePoint(0.16666666666667,   0.78867513459481,   0.03522081090086));
                        break;
                        
                    case 3: //3d, 3 points per dimension
                        mWeights.push_back(0.01497274736603);
                        mWeights.push_back(0.01349962850795);
                        mWeights.push_back(0.00190178826891);
                        mWeights.push_back(0.00760715307442);
                        mWeights.push_back(0.00685871056241);
                        mWeights.push_back(0.00096623512860);
                        mWeights.push_back(0.00024155878219);
                        mWeights.push_back(0.00021779261632);
                        mWeights.push_back(0.00003068198821);
                        mWeights.push_back(0.02395639578565);
                        mWeights.push_back(0.02159940561273);
                        mWeights.push_back(0.00304286123026);
                        mWeights.push_back(0.01217144491907);
                        mWeights.push_back(0.01097393689986);
                        mWeights.push_back(0.00154597620576);
                        mWeights.push_back(0.00038649405150);
                        mWeights.push_back(0.00034846818612);
                        mWeights.push_back(0.00004909118114);
                        mWeights.push_back(0.01497274736603);
                        mWeights.push_back(0.01349962850795);
                        mWeights.push_back(0.00190178826891);
                        mWeights.push_back(0.00760715307442);
                        mWeights.push_back(0.00685871056241);
                        mWeights.push_back(0.00096623512860);
                        mWeights.push_back(0.00024155878219);
                        mWeights.push_back(0.00021779261632);
                        mWeights.push_back(0.00003068198821);
                        
                        mPoints.push_back(MakePoint(0.10000000001607,   0.11270166540000,   0.08872983347426));
                        mPoints.push_back(MakePoint(0.44364916730000,   0.11270166540000,   0.05000000000803));
                        mPoints.push_back(MakePoint(0.78729833458393,   0.11270166540000,   0.01127016654181));
                        mPoints.push_back(MakePoint(0.05635083270000,   0.50000000000000,   0.05000000000803));
                        mPoints.push_back(MakePoint(0.25000000000000,   0.50000000000000,   0.02817541635000));
                        mPoints.push_back(MakePoint(0.44364916730000,   0.50000000000000,   0.00635083269197));
                        mPoints.push_back(MakePoint(0.01270166538393,   0.88729833460000,   0.01127016654181));
                        mPoints.push_back(MakePoint(0.05635083270000,   0.88729833460000,   0.00635083269197));
                        mPoints.push_back(MakePoint(0.10000000001607,   0.88729833460000,   0.00143149884212));
                        mPoints.push_back(MakePoint(0.10000000001607,   0.11270166540000,   0.39364916729197));
                        mPoints.push_back(MakePoint(0.44364916730000,   0.11270166540000,   0.22182458365000));
                        mPoints.push_back(MakePoint(0.78729833458393,   0.11270166540000,   0.05000000000803));
                        mPoints.push_back(MakePoint(0.05635083270000,   0.50000000000000,   0.22182458365000));
                        mPoints.push_back(MakePoint(0.25000000000000,   0.50000000000000,   0.12500000000000));
                        mPoints.push_back(MakePoint(0.44364916730000,   0.50000000000000,   0.02817541635000));
                        mPoints.push_back(MakePoint(0.01270166538393,   0.88729833460000,   0.05000000000803));
                        mPoints.push_back(MakePoint(0.05635083270000,   0.88729833460000,   0.02817541635000));
                        mPoints.push_back(MakePoint(0.10000000001607,   0.88729833460000,   0.00635083269197));
                        mPoints.push_back(MakePoint(0.10000000001607,   0.11270166540000,   0.69856850110968));
                        mPoints.push_back(MakePoint(0.44364916730000,   0.11270166540000,   0.39364916729197));
                        mPoints.push_back(MakePoint(0.78729833458393,   0.11270166540000,   0.08872983347426));
                        mPoints.push_back(MakePoint(0.05635083270000,   0.50000000000000,   0.39364916729197));
                        mPoints.push_back(MakePoint(0.25000000000000,   0.50000000000000,   0.22182458365000));
                        mPoints.push_back(MakePoint(0.44364916730000,   0.50000000000000,   0.05000000000803));
                        mPoints.push_back(MakePoint(0.01270166538393,   0.88729833460000,   0.08872983347426));
                        mPoints.push_back(MakePoint(0.05635083270000,   0.88729833460000,   0.05000000000803));
                        mPoints.push_back(MakePoint(0.10000000001607,   0.88729833460000,   0.01127016654181));
                        break;
                        
                    case 4: //3d, 4 points per dimension
                        mWeights.push_back(0.00423982561968);
                        mWeights.push_back(0.00572288385156);
                        mWeights.push_back(0.00281885467361);
                        mWeights.push_back(0.00031634320391);
                        mWeights.push_back(0.00412036229051);
                        mWeights.push_back(0.00556163317318);
                        mWeights.push_back(0.00273942929295);
                        mWeights.push_back(0.00030742976838);
                        mWeights.push_back(0.00099965677330);
                        mWeights.push_back(0.00134932898618);
                        mWeights.push_back(0.00066462336430);
                        mWeights.push_back(0.00007458670588);
                        mWeights.push_back(0.00002360309872);
                        mWeights.push_back(0.00003185928022);
                        mWeights.push_back(0.00001569255698);
                        mWeights.push_back(0.00000176108183);
                        mWeights.push_back(0.00794866986669);
                        mWeights.push_back(0.01072905315027);
                        mWeights.push_back(0.00528468555374);
                        mWeights.push_back(0.00059306865848);
                        mWeights.push_back(0.00772470439029);
                        mWeights.push_back(0.01042674628127);
                        mWeights.push_back(0.00513578175757);
                        mWeights.push_back(0.00057635807584);
                        mWeights.push_back(0.00187411992466);
                        mWeights.push_back(0.00252967258912);
                        mWeights.push_back(0.00124601155388);
                        mWeights.push_back(0.00013983242583);
                        mWeights.push_back(0.00004425022545);
                        mWeights.push_back(0.00005972861231);
                        mWeights.push_back(0.00002941983138);
                        mWeights.push_back(0.00000330161175);
                        mWeights.push_back(0.00794866986669);
                        mWeights.push_back(0.01072905315027);
                        mWeights.push_back(0.00528468555374);
                        mWeights.push_back(0.00059306865848);
                        mWeights.push_back(0.00772470439029);
                        mWeights.push_back(0.01042674628127);
                        mWeights.push_back(0.00513578175757);
                        mWeights.push_back(0.00057635807584);
                        mWeights.push_back(0.00187411992466);
                        mWeights.push_back(0.00252967258912);
                        mWeights.push_back(0.00124601155388);
                        mWeights.push_back(0.00013983242583);
                        mWeights.push_back(0.00004425022545);
                        mWeights.push_back(0.00005972861231);
                        mWeights.push_back(0.00002941983138);
                        mWeights.push_back(0.00000330161175);
                        mWeights.push_back(0.00423982561968);
                        mWeights.push_back(0.00572288385156);
                        mWeights.push_back(0.00281885467361);
                        mWeights.push_back(0.00031634320391);
                        mWeights.push_back(0.00412036229051);
                        mWeights.push_back(0.00556163317318);
                        mWeights.push_back(0.00273942929295);
                        mWeights.push_back(0.00030742976838);
                        mWeights.push_back(0.00099965677330);
                        mWeights.push_back(0.00134932898618);
                        mWeights.push_back(0.00066462336430);
                        mWeights.push_back(0.00007458670588);
                        mWeights.push_back(0.00002360309872);
                        mWeights.push_back(0.00003185928022);
                        mWeights.push_back(0.00001569255698);
                        mWeights.push_back(0.00000176108183);
                        
                        mPoints.push_back(MakePoint(0.06461106321099,   0.06943184420000,   0.06012499793653));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.06943184420000,   0.04328879995478));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.06943184420000,   0.02132226325621));
                        mPoints.push_back(MakePoint(0.86595709258901,   0.06943184420000,   0.00448606527446));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.33000947820000,   0.04328879995478));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.33000947820000,   0.03116707302848));
                        mPoints.push_back(MakePoint(0.44888729930184,   0.33000947820000,   0.01535160449661));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.33000947820000,   0.00322987757031));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.66999052180000,   0.02132226325621));
                        mPoints.push_back(MakePoint(0.10890625570184,   0.66999052180000,   0.01535160449661));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.66999052180000,   0.00756156217830));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.66999052180000,   0.00159090341870));
                        mPoints.push_back(MakePoint(0.00482078098901,   0.93056815580000,   0.00448606527446));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.93056815580000,   0.00322987757031));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.93056815580000,   0.00159090341870));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.93056815580000,   0.00033471571455));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.06943184420000,   0.28577404826889));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.06943184420000,   0.20575161800155));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.06943184420000,   0.10134469352354));
                        mPoints.push_back(MakePoint(0.86595709258901,   0.06943184420000,   0.02132226325621));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.33000947820000,   0.20575161800155));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.33000947820000,   0.14813706341321));
                        mPoints.push_back(MakePoint(0.44888729930184,   0.33000947820000,   0.07296615908496));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.33000947820000,   0.01535160449661));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.66999052180000,   0.10134469352354));
                        mPoints.push_back(MakePoint(0.10890625570184,   0.66999052180000,   0.07296615908496));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.66999052180000,   0.03594009661688));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.66999052180000,   0.00756156217830));
                        mPoints.push_back(MakePoint(0.00482078098901,   0.93056815580000,   0.02132226325621));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.93056815580000,   0.01535160449661));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.93056815580000,   0.00756156217830));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.93056815580000,   0.00159090341870));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.06943184420000,   0.58018304432012));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.06943184420000,   0.41772022627335));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.06943184420000,   0.20575161800155));
                        mPoints.push_back(MakePoint(0.86595709258901,   0.06943184420000,   0.04328879995478));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.33000947820000,   0.41772022627335));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.33000947820000,   0.30075023588863));
                        mPoints.push_back(MakePoint(0.44888729930184,   0.33000947820000,   0.14813706341321));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.33000947820000,   0.03116707302848));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.66999052180000,   0.20575161800155));
                        mPoints.push_back(MakePoint(0.10890625570184,   0.66999052180000,   0.14813706341321));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.66999052180000,   0.07296615908496));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.66999052180000,   0.01535160449661));
                        mPoints.push_back(MakePoint(0.00482078098901,   0.93056815580000,   0.04328879995478));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.93056815580000,   0.03116707302848));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.93056815580000,   0.01535160449661));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.93056815580000,   0.00322987757031));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.06943184420000,   0.80583209465249));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.06943184420000,   0.58018304432012));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.06943184420000,   0.28577404826889));
                        mPoints.push_back(MakePoint(0.86595709258901,   0.06943184420000,   0.06012499793653));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.33000947820000,   0.58018304432012));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.33000947820000,   0.41772022627335));
                        mPoints.push_back(MakePoint(0.44888729930184,   0.33000947820000,   0.20575161800155));
                        mPoints.push_back(MakePoint(0.62347184427491,   0.33000947820000,   0.04328879995478));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.66999052180000,   0.28577404826889));
                        mPoints.push_back(MakePoint(0.10890625570184,   0.66999052180000,   0.20575161800155));
                        mPoints.push_back(MakePoint(0.22110322249816,   0.66999052180000,   0.10134469352354));
                        mPoints.push_back(MakePoint(0.30709631152509,   0.66999052180000,   0.02132226325621));
                        mPoints.push_back(MakePoint(0.00482078098901,   0.93056815580000,   0.06012499793653));
                        mPoints.push_back(MakePoint(0.02291316667491,   0.93056815580000,   0.04328879995478));
                        mPoints.push_back(MakePoint(0.04651867752509,   0.93056815580000,   0.02132226325621));
                        mPoints.push_back(MakePoint(0.06461106321099,   0.93056815580000,   0.00448606527446));
                        break;
                        
                    default:
                        EXCEPTION("Number of gauss points per dimension not supported.");
                }
            }
            break;
            
            default:
                EXCEPTION("Gauss points not available for this dimension.");
        }
        
    }
    
    /**
     * Get a quadrature point.
     * 
     * @param index The index of the point to return.
     * @return A gaussian quadrature point.
     */
    const c_vector<double, ELEM_DIM>& rGetQuadPoint(int index) const
    {
        assert(index < mNumQuadPoints);
        return mPoints[index];
    }
    
    /**
     * Get the weight associated with a quadrature point.
     */
    double GetWeight(int index) const
    {
        assert(index < mNumQuadPoints);
        return mWeights[index];
    }
    
    /**
     * Get the number of quadrature points. This is the number of points in 
     * each dimension, raised to the power of the number of dimensions.
     */
    int GetNumQuadPoints() const
    {
        return mNumQuadPoints;
    }
    
};

/**
 * Abstract base class for basis functions. There are methods to compute
 * the value and derivative of a particular basis function, or all basis
 * functions on an element together.
 *
 * The methods are documented more fully in the LinearBasisFunction class.
 *
 * @see LinearBasisFunction
 */
template <int ELEM_DIM>
class AbstractBasisFunction
{

public:
    virtual double ComputeBasisFunction(const c_vector<double, ELEM_DIM> &rPoint, int basisIndex) const =0;
    virtual c_vector<double, ELEM_DIM> ComputeBasisFunctionDerivative(const c_vector<double, ELEM_DIM> &rPoint, int basisIndex) const =0;
    virtual c_vector<double, ELEM_DIM+1>& rComputeBasisFunctions(const c_vector<double,ELEM_DIM> &rPoint) =0;
    virtual c_matrix<double, ELEM_DIM, ELEM_DIM+1>& rComputeBasisFunctionDerivatives(const c_vector<double,ELEM_DIM> &rPoint) =0;
    virtual c_matrix<double, ELEM_DIM, ELEM_DIM+1>& rGetTransformedBasisFunctionDerivativesReference(void) =0;
    virtual void UpdateTransformedBasisFunctionDerivatives(const c_vector<double,ELEM_DIM> &rPoint, const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian, bool computeDerivs=true) =0;
    virtual ~AbstractBasisFunction()
    { };
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class AbstractBasisFunction<0>
{
public:
    virtual double ComputeBasisFunction(const Point<0> &rPoint, int basisIndex) const =0;
    virtual c_vector<double, 1> ComputeBasisFunctions(const Point<0> &rPoint) const =0;
    virtual ~AbstractBasisFunction()
    { };
};


template <int ELEM_DIM>
class LinearBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
private:
    c_vector<double, ELEM_DIM+1> mBasisValues;
    c_matrix<double, ELEM_DIM, ELEM_DIM+1> mBasisGradValues;
    c_matrix<double, ELEM_DIM, ELEM_DIM+1> mTransformedBasisGradValues;
public:
    double ComputeBasisFunction(const c_vector<double,ELEM_DIM> &rPoint, int basisIndex) const;
    c_vector<double, ELEM_DIM> ComputeBasisFunctionDerivative(const c_vector<double,ELEM_DIM> &rPoint, int basisIndex) const;
    
    c_vector<double, ELEM_DIM+1>& rComputeBasisFunctions(const c_vector<double,ELEM_DIM> &rPoint);
    c_matrix<double, ELEM_DIM, ELEM_DIM+1>& rComputeBasisFunctionDerivatives(const c_vector<double,ELEM_DIM> &rPoint);
    
    c_matrix<double, ELEM_DIM, ELEM_DIM+1>& rGetTransformedBasisFunctionDerivativesReference();
    void UpdateTransformedBasisFunctionDerivatives(const c_vector<double,ELEM_DIM> &rPoint,
            const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian,
            bool computeDerivs=true);
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class LinearBasisFunction<0> : public AbstractBasisFunction<0>
{
private:
    c_vector<double, 1> mBasisValues;
public:
    double ComputeBasisFunction(const Point<0> &rPoint, int basisIndex) const;
    c_vector<double, 1>& rComputeBasisFunctions(const Point<0> &rPoint);
};



/**
 * Compute a basis function at a point within an element (3d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 *
 * \todo basisIndex should be unsigned (ticket:114)
 */
template <>
double LinearBasisFunction<3>::ComputeBasisFunction(
    const c_vector<double,3> &rPoint,
    int basisIndex) const
{
    assert(basisIndex <= 3);
    assert(basisIndex >= 0);

    switch (basisIndex)
    {
    case 0:
    return 1.0 - rPoint[0] - rPoint[1] - rPoint[2];
    break;
    case 1:
    return rPoint[0];
    break;
    case 2:
    return rPoint[1];
    break;
    case 3:
    return rPoint[2];
    break;
    default:
    ; //not possible to get here because of assertions above
    }
    return 0.0; // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element (2d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 *
 * \todo basisIndex should be unsigned (ticket:114)
 */
template <>
double LinearBasisFunction<2>::ComputeBasisFunction(
    const c_vector<double,2> &rPoint,
    int basisIndex) const
{
    assert(basisIndex <= 2);
    assert(basisIndex >= 0);

    switch (basisIndex)
    {
    case 0:
    return 1.0 - rPoint[0] - rPoint[1];
    break;
    case 1:
    return rPoint[0];
    break;
    case 2:
    return rPoint[1];
    break;
    default:
    ; //not possible to get here because of assertions above
    }
    return 0.0; // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element (1d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 *
 * \todo basisIndex should be unsigned (ticket:114)
 */
template <>
double LinearBasisFunction<1>::ComputeBasisFunction(
    const c_vector<double,1> &rPoint,
    int basisIndex) const
{
    assert(basisIndex <= 1);
    assert(basisIndex >= 0);

    switch (basisIndex)
    {
    case 0:
    return 1.0 - rPoint[0];
    break;
    case 1:
    return rPoint[0];
    break;
    default:
    ; //not possible to get here because of assertions above
    }
    return 0.0; // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element (0d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
double LinearBasisFunction<0>::ComputeBasisFunction(const Point<0> &rPoint, int basisIndex) const
{
    assert(basisIndex == 0);
    return 1.0;
}

/**
 * Compute the derivative of a basis function at a point within a
 * canonical element (3d case).
 *
 * @param rPoint (unused) The point at which to compute the basis function.
 *     The results are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector
 *     (c_vector<double, ELEM_DIM> instance) giving the derivative
 *     along each axis.
 *
 * \todo basisIndex should be unsigned (ticket:114)
 */
template <>
c_vector<double, 3> LinearBasisFunction<3>::ComputeBasisFunctionDerivative(
    const c_vector<double,3>&,
    int basisIndex) const
{
    assert(basisIndex <= 3);
    assert(basisIndex >= 0);

    c_vector<double, 3> gradN;
    switch (basisIndex)
    {
    case 0:
    gradN(0) = -1;
    gradN(1) = -1;
    gradN(2) = -1;
    break;
    case 1:
    gradN(0) =  1;
    gradN(1) =  0;
    gradN(2) =  0;
    break;
    case 2:
    gradN(0) =  0;
    gradN(1) =  1;
    gradN(2) =  0;
    break;
    case 3:
    gradN(0) =  0;
    gradN(1) =  0;
    gradN(2) =  1;
    break;
    default:
    ; //not possible to get here because of assertions above
    }
    return gradN;
}

/**
 * Compute the derivative of a basis function at a point within a
 * canonical element (2d case).
 *
 * @param rPoint (unused) The point at which to compute the basis function.
 *     The results are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector
 *     (c_vector<double, ELEM_DIM> instance) giving the derivative
 *     along each axis.
 *
 * \todo basisIndex should be unsigned (ticket:114)
 */
template <>
c_vector<double, 2> LinearBasisFunction<2>::ComputeBasisFunctionDerivative(
    const c_vector<double,2>&,
    int basisIndex) const
{
    assert(basisIndex <= 2);
    assert(basisIndex >= 0);

    c_vector<double, 2> gradN;
    switch (basisIndex)
    {
    case 0:
    gradN(0) = -1;
    gradN(1) = -1;
    break;
    case 1:
    gradN(0) =  1;
    gradN(1) =  0;
    break;
    case 2:
    gradN(0) =  0;
    gradN(1) =  1;
    break;
    default:
    ; //not possible to get here because of assertions above
    }
    return gradN;
}

/**
 * Compute the derivative of a basis function at a point within a
 * canonical element (1d case).
 *
 * @param rPoint (unused) The point at which to compute the basis function.
 *     The results are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector
 *     (c_vector<double, ELEM_DIM> instance) giving the derivative
 *     along each axis.
 *
 * \todo basisIndex should be unsigned (ticket:114)
 */
template <>
c_vector<double, 1> LinearBasisFunction<1>::ComputeBasisFunctionDerivative(
    const c_vector<double,1>&,
    int basisIndex) const
{
    assert(basisIndex <= 1);
    assert(basisIndex >= 0);

    c_vector<double, 1> gradN;
    switch (basisIndex)
    {
    case 0:
    gradN(0) = -1;
    break;
    case 1:
    gradN(0) =  1;
    break;
    default:
    ; //not possible to get here because of assertions above
    }
    return gradN;
}



/**
 * Compute all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @return The values of the basis functions, in local index order.
 */
template <int ELEM_DIM>
c_vector<double, ELEM_DIM+1>& LinearBasisFunction<ELEM_DIM>::rComputeBasisFunctions(const c_vector<double,ELEM_DIM> &rPoint)
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
//    c_vector<double, ELEM_DIM+1> basisValues;
    for (int i=0; i<ELEM_DIM+1; i++)
    {
        mBasisValues(i) = ComputeBasisFunction(rPoint, i);
    }
    return mBasisValues;
}

/**
 * Compute all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @return The values of the basis functions, in local index order.
 *
 */
c_vector<double, 1>& LinearBasisFunction<0>::rComputeBasisFunctions(const Point<0> &rPoint)
{
//    c_vector<double, 1> basisValues;
    mBasisValues(0) = ComputeBasisFunction(rPoint, 0);
    return mBasisValues;
}

/**
 * Compute the derivatives of all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @return The derivatives of the basis functions as the column vectors of
 *     a matrix in local index order.
 */
template <int ELEM_DIM>
c_matrix<double, ELEM_DIM, ELEM_DIM+1>& LinearBasisFunction<ELEM_DIM>::rComputeBasisFunctionDerivatives(const c_vector<double,ELEM_DIM> &rPoint)
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
//    c_matrix<double, ELEM_DIM, ELEM_DIM+1> basisGradValues;
    
    for (unsigned j=0;j<ELEM_DIM+1;j++)
    {
        matrix_column<c_matrix<double, ELEM_DIM, ELEM_DIM+1> > column(mBasisGradValues, j);
        column = ComputeBasisFunctionDerivative(rPoint, j);
    }
    
    return mBasisGradValues;
}


/**
 * Compute the derivatives of all basis functions at a point within an element.
 * This method will transform the results, for use within gaussian quadrature
 * for example.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @param inverseJacobian The inverse of the Jacobian matrix mapping the real
 *     element into the canonical element.
 * @return The derivatives of the basis functions, in local index order. Each
 *     entry is a vector (c_vector<double, SPACE_DIM> instance) giving the
 *     derivative along each axis.
 */
template <int ELEM_DIM>
void LinearBasisFunction<ELEM_DIM>::UpdateTransformedBasisFunctionDerivatives(const c_vector<double,ELEM_DIM> &rPoint, const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian, bool computeDerivs)
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
//    c_matrix<double, ELEM_DIM, ELEM_DIM+1> basisGradValues = ComputeBasisFunctionDerivatives(rPoint);
    if (computeDerivs)
    {
        rComputeBasisFunctionDerivatives(rPoint);
    }
    mTransformedBasisGradValues = prod(trans(rInverseJacobian), mBasisGradValues);
}

template <int ELEM_DIM>
c_matrix<double, ELEM_DIM, ELEM_DIM+1>& LinearBasisFunction<ELEM_DIM>::rGetTransformedBasisFunctionDerivativesReference(void)
{
    return mTransformedBasisGradValues;
}


/**
 *  BidomainDg0Assembler
 *
 *  inherits from AbstractLinearDynamicProblemAssembler<ELEM_DIM, SPACE_DIM, 2> (the
 *  2 representing the number of unknowns (ie voltage and extracellular potential)).
 *
 *  This assembler interpolates quantities such as ionic currents and stimuli from
 *  their nodal values (obtained from a BidomainPde) onto a gauss point, and uses
 *  the interpolated values in assembly. The assembler also creates boundary conditions,
 *  which are zero-Neumann boundary conditions on the surface unless
 *  SetFixedExtracellularPotentialNodes() is called.
 *
 *  The user should call Solve() from the superclass AbstractLinearDynamicProblemAssembler.
 *
 *  NOTE: if any cells have a non-zero extracellular stimulus, phi_e must be fixed at some
 *  nodes (using SetFixedExtracellularPotentialNodes() ), else no solution is possible.
 */

template<int ELEMENT_DIM, int SPACE_DIM>
class BidomainDg0Assembler //: public AbstractLinearDynamicProblemAssembler<ELEMENT_DIM, SPACE_DIM, 2>
{
protected :
    //
    // From AbstractAssembler
    //
    bool mWeAllocatedBasisFunctionMemory;
    
    /*< Mesh to be solved on */
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;
    
    /*< Boundary conditions to be applied */
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditions;
    
    /*< Basis function for use with normal elements */
    AbstractBasisFunction<ELEMENT_DIM> *mpBasisFunction;
    /*< Basis function for use with boundary elements */
    AbstractBasisFunction<ELEMENT_DIM-1> *mpSurfaceBasisFunction;
    
    /*< Quadrature rule for use on normal elements */
    GaussianQuadratureRule<ELEMENT_DIM> *mpQuadRule;
    /*< Quadrature rule for use on boundary elements */
    GaussianQuadratureRule<ELEMENT_DIM-1> *mpSurfaceQuadRule;

    /**
     *  The CURRENT SOLUTION as a replicated vector for linear dynamic problems. 
     *  (Empty for a static problem). The CURRENT GUESS for nonlinear problems
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;
    

    /*< bool stating whether the problem is a linear or nonlinear one */
    bool mProblemIsLinear; 

    /** 
     *  The linear system that is assembled in linear pde problems. Not used in
     *  nonlinear problems
     */
    LinearSystem *mpLinearSystem;
    
    /**
     * mMatrixIsConstant is a flag to say whether the matrix of the system
     * needs to be assembled at each time step. (Linear problems only).
     */
    bool mMatrixIsConstant;
    
    /**
     * mMatrixIsAssembled is a flag to say whether the matrix has been assembled 
     * for the current time step. (Linear problems only).
     */
    bool mMatrixIsAssembled;
    
    //
    // From AbstractLinearAssembler
    //
    /**
     * The linear solver used to solve the linear system at each time step.
     */
    AbstractLinearSolver *mpLinearSolver;
    bool mWeAllocatedSolverMemory;

    //
    // From AbstractLinearDynamicProblemAssembler
    //
    double mTstart;
    double mTend;
    double mDt, mDtInverse;
    
    bool   mTimesSet;
    bool   mInitialConditionSet;
    
    Vec    mInitialCondition;
    

private:
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    
    // quantities to be interpolated
    double mIionic;
    double mIIntracellularStimulus;
    double mIExtracellularStimulus;
    
    bool mNullSpaceCreated;
    
    std::vector<unsigned> mFixedExtracellularPotentialNodes;
    
    
    void ResetInterpolatedQuantities( void )
    {
        mIionic=0;
        mIIntracellularStimulus=0;
        mIExtracellularStimulus=0;
    }
    
    
    void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM>* pNode)
    {
        unsigned node_global_index = pNode->GetIndex();
        
        mIionic                 += phi_i * mpBidomainPde->GetIionicCacheReplicated()[ node_global_index ];
        mIIntracellularStimulus += phi_i * mpBidomainPde->GetIntracellularStimulusCacheReplicated()[ node_global_index ];
        mIExtracellularStimulus += phi_i * mpBidomainPde->GetExtracellularStimulusCacheReplicated()[ node_global_index ];
    }
    
    /** 
     *  ComputeMatrixTerm()
     * 
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness matrix.
     */        
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        c_vector<double, SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */)
    {
        // get bidomain parameters
        double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
        double Cm = mpBidomainPde->GetCapacitance();
        
        c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_i = mpBidomainPde->GetIntracellularConductivityTensor();
        c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_e = mpBidomainPde->GetExtracellularConductivityTensor();
        
        
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
            prod(trans(rGradPhi), temp);
            
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
            outer_prod(rPhi, rPhi);
            
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
            prod(trans(rGradPhi), temp2);
            
            
        c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret;
        
        // even rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice00(ret, slice (0, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        slice00 = (Am*Cm/this->mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi ;
        
        // odd rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice10(ret, slice (1, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        slice10 = grad_phi_sigma_i_grad_phi;
        
        // even rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice01(ret, slice (0, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice01 = grad_phi_sigma_i_grad_phi;
        
        // odd rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice11(ret, slice (1, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice11 = grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi;
        
        return ret;
    }
    
    
    /**
     *  ComputeVectorTerm()
     * 
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness vector.
     */
    virtual c_vector<double,2*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        c_vector<double,SPACE_DIM>  &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */)
    {
        // get bidomain parameters
        double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
        double Cm = mpBidomainPde->GetCapacitance();
        
        c_vector<double,2*(ELEMENT_DIM+1)> ret;
        
        vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_V  (ret, slice (0, 2, ELEMENT_DIM+1));
        vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEMENT_DIM+1)); 
        
        // u(0) = voltage
        slice_V   =  (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus) * rPhi;
        slice_Phi =  -mIExtracellularStimulus * rPhi;
        
        return ret;
    }
    
    
    
    
    /** ComputeSurfaceLhsTerm()
     * 
     *  This method is called by AssembleOnSurfaceElement() and tells the 
     *  assembler what to add to the element stiffness matrix arising 
     *  from surface element contributions.
     * 
     *  NOTE: this method has to be implemented but shouldn't ever be called -
     *  because all bidomain problems (currently) just have zero Neumann boundary
     *  conditions and the AbstractLinearAssmebler::AssembleSystem() method
     *  will realise this and not loop over surface elements.
     */
#define COVERAGE_IGNORE //see NOTE above
    virtual c_vector<double, 2*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double,ELEMENT_DIM> &rPhi,
        c_vector<double,SPACE_DIM> &rX)
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_grad_v_dot_n     = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, Point<SPACE_DIM>(rX), 0);
        double D_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, Point<SPACE_DIM>(rX), 1);
        
        c_vector<double, 2*ELEMENT_DIM> ret;
        for (int i=0; i<ELEMENT_DIM; i++)
        {
            ret(2*i)   = rPhi(i)*D_times_grad_v_dot_n;
            ret(2*i+1) = rPhi(i)*D_times_grad_phi_e_dot_n;
        }
        
        return ret;
    }
#undef COVERAGE_IGNORE
    
    
    
    
    /**
     *  PrepareForAssembleSystem
     * 
     *  Called at the beginning of AbstractLinearAssmebler::AssembleSystem() 
     *  after the system. Here, used to integrate cell odes.
     */
    virtual void PrepareForAssembleSystem(Vec currentSolution, double time)
    {
        mpBidomainPde->SolveCellSystems(currentSolution, time, time+this->mDt);
    }
    
    /**
     *  FinaliseAssembleSystem
     * 
     *  Called by AbstractLinearAssmebler::AssembleSystem() after the system
     *  has been assembled. Here, used to set up a null basis.
     */
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime)
    {
        // if there are no fixed nodes then set up the null basis.
        if ( (mFixedExtracellularPotentialNodes.size()==0) && (!mNullSpaceCreated) )
        {
            //create null space for matrix and pass to linear system
            Vec nullbasis[1];
            unsigned lo, hi;
            
            mpBidomainPde->GetOwnershipRange(lo, hi);
            VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*this->mpMesh->GetNumNodes(), &nullbasis[0]);
            double* p_nullbasis;
            VecGetArray(nullbasis[0], &p_nullbasis);
            
            for (unsigned global_index=lo; global_index<hi; global_index++)
            {
                unsigned local_index = global_index - lo;
                p_nullbasis[2*local_index  ] = 0;
                p_nullbasis[2*local_index+1] = 1;
            }
            VecRestoreArray(nullbasis[0], &p_nullbasis);
            VecAssemblyBegin(nullbasis[0]);
            VecAssemblyEnd(nullbasis[0]);
            
            this->mpLinearSystem->SetNullBasis(nullbasis, 1);
            
            VecDestroy(nullbasis[0]);
            mNullSpaceCreated = true;
        }
    }
    
    
public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     */
    BidomainDg0Assembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         BidomainPde<SPACE_DIM>* pPde,
                         int numQuadPoints = 2,
                         double linearSolverRelativeTolerance = 1e-6)
    {
        std::cout << "In pointless constructor." << std::endl;
        //
        // From AbstractAssembler
        //
        
        // Initialise mesh and bcs to null, so we can check they
        // have been set before attempting to solve
        mpMesh = NULL;
        mpBoundaryConditions = NULL;
        
        mWeAllocatedBasisFunctionMemory = false; // sic
        LinearBasisFunction<ELEMENT_DIM> *pBasisFunction = new LinearBasisFunction<ELEMENT_DIM>();
        LinearBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction = new LinearBasisFunction<ELEMENT_DIM-1>();
        SetBasisFunctions(pBasisFunction, pSurfaceBasisFunction);
        mWeAllocatedBasisFunctionMemory = true;
        
        mpQuadRule = NULL;
        mpSurfaceQuadRule = NULL;
        SetNumberOfQuadraturePointsPerDimension(numQuadPoints);

        mMatrixIsAssembled = false;
        
        //
        // From AbstractLinearAssembler
        //
        mpLinearSolver = new SimpleLinearSolver(linearSolverRelativeTolerance);
        mWeAllocatedSolverMemory = true;
        
        this->mpLinearSystem = NULL;
        this->mMatrixIsConstant = false;
        this->mMatrixIsAssembled = false;
        
        this->mProblemIsLinear = true;
        
        //
        // From us
        //
        assert(pPde != NULL);
        assert(pMesh != NULL);
        
        mpBidomainPde = pPde;
        this->mpMesh = pMesh;
        
        // set up boundary conditions
        this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>( this->mpMesh->GetNumNodes() );
        
        // define zero neumann boundary conditions everywhere
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,0); // first unknown, ie voltage
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,1); // second unknown, ie phi_e
        
        this->mMatrixIsAssembled = false;
        mNullSpaceCreated = false;
        
        this->SetMatrixIsConstant();
        
        mFixedExtracellularPotentialNodes.resize(0);
    }
    
    
    virtual ~BidomainDg0Assembler()
    {
        //
        // From AbstractAssembler
        //
        
        // Basis functions, if we used the default.
        if (mWeAllocatedBasisFunctionMemory)
        {
            delete mpBasisFunction;
            delete mpSurfaceBasisFunction;
            mWeAllocatedBasisFunctionMemory = false;
        }
        
        // Quadrature rules
        if (mpQuadRule) delete mpQuadRule;
        if (mpSurfaceQuadRule) delete mpSurfaceQuadRule;
        
        //
        // From AbstractLinearAssembler
        //
        if (this->mpLinearSystem != NULL)
        {
            delete this->mpLinearSystem;
        }
        
        this->mpLinearSystem=NULL;
        
        if(mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
        
        // From us
        delete this->mpBoundaryConditions;
    }
    
    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to 
     *  zero. This does not necessarily have to be called. If it is not, phi_e 
     *  is only defined up to a constant.
     * 
     *  @param the nodes to be fixed.
     * 
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes)
    {
        assert(fixedExtracellularPotentialNodes.size() > 0);
        for (unsigned i=0; i<fixedExtracellularPotentialNodes.size(); i++)
        {
            if ( (int) fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
            {
                EXCEPTION("Fixed node number must be less than total number nodes");
            }
        }
        
        mFixedExtracellularPotentialNodes = fixedExtracellularPotentialNodes;
        
        for (unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition
            = new ConstBoundaryCondition<SPACE_DIM>(0.0);
            
            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(mFixedExtracellularPotentialNodes[i]);
            
            this->mpBoundaryConditions->AddDirichletBoundaryCondition(p_node, p_boundary_condition, 1);
        }
    }
    
    //
    // From AbstractLinearDynamicProblemAssembler
    //
    /**
     *  Set the times to solve between, and the time step to use
     */
    void SetTimes(double Tstart, double Tend, double dt)
    {
        mTstart = Tstart;
        mTend   = Tend;
        mDt     = dt;
        mDtInverse = 1/dt;
        
        if (mTstart >= mTend)
        {
            EXCEPTION("Starting time has to less than ending time");
        }
        if (mDt <= 0)
        {
            EXCEPTION("Time step has to be greater than zero");
        }
        
        assert(mDt <= mTend - mTstart + 1e-10);
        
        mTimesSet = true;
    }
    
    /**
     *  Set the initial condition
     */
    void SetInitialCondition(Vec initCondition)
    {
        mInitialCondition = initCondition;
        mInitialConditionSet = true;
    }
    
    
    /**
     *  Solve a dynamic PDE over the time period specified through SetTimes()
     *  and the initial conditions specified through SetInitialCondition().
     * 
     *  SetTimes() and SetInitialCondition() must be called before Solve(), and 
     *  the mesh and pde must have been set.
     */
    Vec Solve()
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        this->PrepareForSolve();
        
        double t = mTstart;
        Vec currentSolution = mInitialCondition;
        Vec nextSolution;
        while ( t < mTend - 1e-10 )
        {
            this->AssembleSystem(currentSolution, t);
            
            nextSolution = this->mpLinearSystem->Solve(this->mpLinearSolver);
            
            t += mDt;
            // Avoid memory leaks
            if (currentSolution != mInitialCondition)
            {
                VecDestroy(currentSolution);
            }
            currentSolution = nextSolution;
            
        }
        
        return currentSolution;
    }
    
    
    //
    // From AbstractLinearSolver
    //
    /**
     *  Manually re-set the linear system solver (which by default 
     *  is a SimpleLinearSolver)
     */
    void SetLinearSolver(AbstractLinearSolver *pLinearSolver)
    {
        if(mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
        mpLinearSolver = pLinearSolver;
        
        // make sure new solver knows matrix is constant
        if (this->mMatrixIsConstant)
        {
            SetMatrixIsConstant();
        }
    }
    
    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
     */
    void SetMatrixIsConstant()
    {
        this->mMatrixIsConstant = true;
        this->mpLinearSolver->SetMatrixIsConstant();
    }
    
    //
    // From AbstractAssembler
    //
        /**
     *  Calculate the contribution of a single element to the linear system.
     * 
     *  @param rElement The element to assemble on.
     *  @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     *  @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     *  @param currentSolutionOrGuess For the parabolic linear case, the solution at the current 
     *     timestep. NULL for the static linear case. In the nonlinear case, the current
     *     guess.
     *  @param assembleVector a bool stating whether to assemble the load vector (in the 
     *     linear case) or the residual vector (in the nonlinear case)
     *  @param assembleMatrix a bool stating whether to assemble the stiffness matrix (in 
     *     the linear case) or the Jacobian matrix (in the nonlinear case)
     * 
     *  Called by AssembleSystem()
     *  Calls ComputeMatrixTerm() etc
     */
    virtual void AssembleOnElement( Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) > &rAElem,
                                    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> &rBElem,
                                    bool assembleVector,
                                    bool assembleMatrix)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(this->mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(this->mpBasisFunction);
            
            
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function.
         */
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *p_inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (! (mProblemIsLinear && mMatrixIsAssembled) ) // don't need to construct grad_phi or grad_u in that case
        {
            p_inverse_jacobian = rElement.GetInverseJacobian();
            rAElem.clear();
        }
        
        rBElem.clear();
        
        
        const int num_nodes = rElement.GetNumNodes();
        
        // loop over Gauss points
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const c_vector<double, ELEMENT_DIM>& quad_point = quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1>& phi = rBasisFunction.rComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& grad_phi = rBasisFunction.rGetTransformedBasisFunctionDerivativesReference();

            if (! (mProblemIsLinear && mMatrixIsAssembled) ) // don't need to construct grad_phi or grad_u in that case
            {
                rBasisFunction.UpdateTransformedBasisFunctionDerivatives(quad_point, *p_inverse_jacobian, false);
            }
            
            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            c_vector<double,SPACE_DIM> x;
            x.clear();
            
            c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
            c_matrix<double,PROBLEM_DIM,SPACE_DIM> grad_u = zero_matrix<double>(PROBLEM_DIM,SPACE_DIM);
            
            // allow the concrete version of the assembler to interpolate any
            // desired quantities
            ResetInterpolatedQuantities();
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            for (int i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
                const c_vector<double, SPACE_DIM> node_loc = p_node->rGetLocation();
                
                // interpolate x
                x += phi(i)*node_loc;
                
                // interpolate u and grad u if a current solution or guess exists
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                if (mCurrentSolutionOrGuessReplicated.size()>0)
                {
                    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
                    {
                        // If we have a current solution (e.g. this is a dynamic problem)
                        // get the value in a usable form.
                        
                        // NOTE - currentSolutionOrGuess input is actually now redundant at this point -
                        
                        // NOTE - following assumes that, if say there are two unknowns u and v, they
                        // are stored in the curren solution vector as
                        // [U1 V1 U2 V2 ... U_n V_n]
                        u(index_of_unknown) += phi(i)*this->mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*node_global_index + index_of_unknown];

                        if (! (mProblemIsLinear && mMatrixIsAssembled) ) // don't need to construct grad_phi or grad_u in that case
                        {
                            for(unsigned j=0; j<SPACE_DIM; j++)
                            {
                               grad_u(index_of_unknown,j) += grad_phi(j,i)*this->mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*node_global_index + index_of_unknown];
                            }
                        }
                    }        
                }
                
                // allow the concrete version of the assembler to interpolate any
                // desired quantities
                IncrementInterpolatedQuantities(phi(i), p_node);
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            
            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            if(assembleMatrix) 
            {
                noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
            
            if(assembleVector)
            {
                noalias(rBElem) += ComputeVectorTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
        }
    }
    
    
    
    /**
     * Calculate the contribution of a single surface element with Neumann
     * boundary condition to the linear system.
     * 
     * @param rSurfaceElement The element to assemble on.
     * @param rBSurfElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                          c_vector<double, PROBLEM_DIM*ELEMENT_DIM> &rBSurfElem)
    {
        GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
            *(this->mpSurfaceQuadRule);
        AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
            *(this->mpSurfaceBasisFunction);
            
        double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
        
        rBSurfElem.clear();
        
        // loop over Gauss points
        for (int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const c_vector<double, ELEMENT_DIM-1>& quad_point=quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM>& phi = rBasisFunction.rComputeBasisFunctions(quad_point);
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            
            // Location of the gauss point in the original element will be
            // stored in x
            c_vector<double,SPACE_DIM> x;
            x.clear();
                       
            ResetInterpolatedQuantities();
            for (int i=0; i<rSurfaceElement.GetNumNodes(); i++)
            {
                const c_vector<double, SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetLocation();
                x += phi(i)*node_loc;
                
                // allow the concrete version of the assembler to interpolate any
                // desired quantities
                IncrementInterpolatedQuantities(phi(i), rSurfaceElement.GetNode(i));
                
                ///\todo: add interpolation of u as well
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            
            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            ///\todo Improve efficiency of Neumann BC implementation.
            noalias(rBSurfElem) += ComputeVectorSurfaceTerm(rSurfaceElement, phi, x) * wJ;
        }
    }
    
    /**
     *  AssembleSystem - the major method for all assemblers
     * 
     *  Assemble the linear system for a linear PDE, or the residual or Jacobian for
     *  nonlinear PDEs. Loops over each element (and each each surface element if 
     *  there are non-zero Neumann boundary conditions and calls AssembleOnElement() 
     *  and adds the contribution to the linear system.
     * 
     *  Takes in current solution and time if necessary but only used if the problem 
     *  is a dynamic one. This method uses PROBLEM_DIM and can assemble linear systems 
     *  for any number of unknown variables.
     * 
     *  Called by Solve()
     *  Calls AssembleOnElement()
     * 
     *  @param currentSolutionOrGuess The current solution in a linear dynamic problem, 
     *     or the current guess in a nonlinear problem. Should be NULL for linear static 
     *     problems.
     * 
     *  @param currentTime The current time for dynamic problems. Not used in static 
     *     problems.
     * 
     *  @param residualVector The residual vector to be assembled in nonlinear problems
     *     (eg created by the Petsc nonlinear solver). Should be NULL for linear problems.
     * 
     *  @param pJacobianMatrix (A pointer to) the Jacobian matrix to be assembled in 
     *     nonlinear problems (eg created by the Petsc nonlinear solver). Should be 
     *     NULL for linear problems.
     */
    virtual void AssembleSystem(Vec currentSolutionOrGuess=NULL, double currentTime=0.0, Vec residualVector=NULL, Mat* pJacobian=NULL)
    {
        // if a linear problem there mustn't be a residual or jacobian specified
        // otherwise one of them MUST be specifed
        assert(    (mProblemIsLinear && !residualVector && !pJacobian) 
                || (!mProblemIsLinear && (residualVector || pJacobian) ) );
        
        // if the problem is nonlinear the currentSolutionOrGuess MUST be specifed
        assert( mProblemIsLinear || (!mProblemIsLinear && currentSolutionOrGuess ) );
                        
        // Replicate the current solution and store so can be used in
        // AssembleOnElement
        if (currentSolutionOrGuess != NULL)
        {
            this->mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolutionOrGuess);
        }
        
        // the AssembleOnElement type methods will determine if a current solution or
        // current guess exists by looking at the size of the replicated vector, so 
        // check the size is zero if there isn't a current solution
        assert(    ( currentSolutionOrGuess && mCurrentSolutionOrGuessReplicated.size()>0)
                || ( !currentSolutionOrGuess && mCurrentSolutionOrGuessReplicated.size()==0));
        

        // the concrete class can override this following method if there is
        // work to be done before assembly
        PrepareForAssembleSystem(currentSolutionOrGuess, currentTime);
        
        int lo, hi;
        
        if(mProblemIsLinear)
        {
            // linear problem - set up the Linear System if necessary, otherwise zero
            // it.
            if (mpLinearSystem == NULL)
            {
                if (currentSolutionOrGuess == NULL)
                {
                    // static problem, create linear system using the size
                    unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
                    mpLinearSystem = new LinearSystem(size);
                }
                else
                {
                    // use the currrent solution (ie the initial solution)
                    // as the template in the alternative constructor of
                    // LinearSystem. This appears to avoid problems with
                    // VecScatter.
                    mpLinearSystem = new LinearSystem(currentSolutionOrGuess);
                }
            }
            else
            {
                if (mMatrixIsConstant && mMatrixIsAssembled)
                {
                    mpLinearSystem->ZeroRhsVector();
                }
                else
                {
                    mpLinearSystem->ZeroLinearSystem();
                    mMatrixIsAssembled = false;
                }
            }
        }
        else
        {   
            // nonlinear problem - zero residual or jacobian depending on which has
            // been asked for     
            if(residualVector)
            {
                int size;
                VecGetSize(residualVector,&size);
                assert(size==PROBLEM_DIM * this->mpMesh->GetNumNodes());
            
                // Set residual vector to zero
                PetscScalar zero = 0.0;
#if (PETSC_VERSION_MINOR == 2) //Old API
                PETSCEXCEPT( VecSet(&zero, residualVector) );
#else
                PETSCEXCEPT( VecSet(residualVector, zero) );
#endif
            }
            else 
            {
                int size1, size2;
                MatGetSize(*pJacobian,&size1,&size2);
                assert(size1==PROBLEM_DIM * this->mpMesh->GetNumNodes());
                assert(size2==PROBLEM_DIM * this->mpMesh->GetNumNodes());
   
                // Set all entries of jacobian to 0
                MatZeroEntries(*pJacobian);
            }        
        
            // Get our ownership range
            VecGetOwnershipRange(currentSolutionOrGuess, &lo, &hi);
        }
        
                 
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
            iter = this->mpMesh->GetElementIteratorBegin();
        
        // Assume all elements have the same number of nodes...
        const int num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> b_elem;
        

        // decide what we want to assemble. 
        bool assemble_vector = ((mProblemIsLinear) || ((!mProblemIsLinear) && (residualVector!=NULL)));
        bool assemble_matrix = ( (mProblemIsLinear && !mMatrixIsAssembled) || ((!mProblemIsLinear) && (pJacobian!=NULL)) );
       
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            AssembleOnElement(element, a_elem, b_elem, assemble_vector, assemble_matrix);
            
            for (int i=0; i<num_elem_nodes; i++)
            {
                int node1 = element.GetNodeGlobalIndex(i);
                                
                if (assemble_matrix)
                {                    
                    for (int j=0; j<num_elem_nodes; j++)
                    {
                        int node2 = element.GetNodeGlobalIndex(j);
                        
                        for (int k=0; k<PROBLEM_DIM; k++)
                        {
                            for (int m=0; m<PROBLEM_DIM; m++)
                            {
                                if(mProblemIsLinear)
                                {  
                                    // the following expands to, for (eg) the case of two unknowns:
                                    // mpLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                                    // mpLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                                    // mpLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                                    // mpLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                                    mpLinearSystem->AddToMatrixElement( PROBLEM_DIM*node1+k,
                                                                        PROBLEM_DIM*node2+m,
                                                                        a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+m) );
                                }
                                else 
                                {
                                    assert(pJacobian!=NULL); // extra check
                                           
                                    int matrix_index_1 = PROBLEM_DIM*node1+k;
                                    if (lo<=matrix_index_1 && matrix_index_1<hi)
                                    {
                                        int matrix_index_2 = PROBLEM_DIM*node2+m;
                                        PetscScalar value = a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+m);
                                        MatSetValue(*pJacobian, matrix_index_1, matrix_index_2, value, ADD_VALUES);                                
                                    }
                                }
                            }
                        }
                    }
                }

                if(assemble_vector)
                {
                    for (int k=0; k<PROBLEM_DIM; k++)
                    {
                        if(mProblemIsLinear)
                        {
                            mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1+k,b_elem(PROBLEM_DIM*i+k));
                        }
                        else 
                        {
                            assert(residualVector!=NULL); // extra check

                            int matrix_index = PROBLEM_DIM*node1+k;
                            //Make sure it's only done once
                            if (lo<=matrix_index && matrix_index<hi)
                            {
                                PetscScalar value = b_elem(PROBLEM_DIM*i+k);
                                PETSCEXCEPT( VecSetValue(residualVector,matrix_index,value,ADD_VALUES) );
                            }
                        }
                    }
                }
            }
            iter++;
        }
                
        // add the integrals associated with Neumann boundary conditions to the linear system
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator
        surf_iter = this->mpMesh->GetBoundaryElementIteratorBegin();
        

        ////////////////////////////////////////////////////////
        // loop over surface elements
        ////////////////////////////////////////////////////////

        // note, the following condition is not true of Bidomain or Monodomain
        if (this->mpBoundaryConditions->AnyNonZeroNeumannConditions()==true)
        {
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const int num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
                    
                    ///\todo Check surf_element is in the Neumann surface in an efficient manner
                    /// e.g. by iterating over boundary conditions!
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        AssembleOnSurfaceElement(surf_element, b_surf_elem);
                        
                        for (int i=0; i<num_surf_nodes; i++)
                        {
                            int node_index = surf_element.GetNodeGlobalIndex(i);
                            
                            for (int k=0; k<PROBLEM_DIM; k++)
                            {
                                if(mProblemIsLinear)
                                {
                                    mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node_index + k, b_surf_elem(PROBLEM_DIM*i+k));
                                }
                                else if(residualVector!=NULL)
                                {
                                    int matrix_index = PROBLEM_DIM*node_index + k;

                                    PetscScalar value = b_surf_elem(PROBLEM_DIM*i+k);
                                    if (lo<=matrix_index && matrix_index<hi)
                                    {
                                        PETSCEXCEPT( VecSetValue(residualVector, matrix_index, value, ADD_VALUES) );
                                    }
                                }
                            }
                        }
                    }
                    surf_iter++;
                }
            }
        }

        
        if(mProblemIsLinear)
        {
            if (mMatrixIsAssembled)
            {
                mpLinearSystem->AssembleRhsVector();
            }
            else
            {
                mpLinearSystem->AssembleIntermediateLinearSystem();
            }
        }
        else if(pJacobian)
        {
            MatAssemblyBegin(*pJacobian, MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(*pJacobian, MAT_FLUSH_ASSEMBLY);
        }
        
        
        // Apply dirichlet boundary conditions
        if(mProblemIsLinear)
        {
            this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*mpLinearSystem, mMatrixIsAssembled);
        }
        else if(residualVector)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentSolutionOrGuess, residualVector);
        }        
        else if(pJacobian)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pJacobian);
        }
        
        
                    
        if(mProblemIsLinear)
        {        
            if (mMatrixIsAssembled)
            {
                mpLinearSystem->AssembleRhsVector();
            }
            else
            {
                mpLinearSystem->AssembleFinalLinearSystem();
            }
            mMatrixIsAssembled = true;
        }
        else if(residualVector)
        {
            VecAssemblyBegin(residualVector);
            VecAssemblyEnd(residualVector);
        }
        else if(pJacobian)
        {
            MatAssemblyBegin(*pJacobian, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(*pJacobian, MAT_FINAL_ASSEMBLY);
        }        
        
        // overload this method if the assembler has to do anything else
        // required (like setting up a null basis (see BidomainDg0Assembler))
        FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
    }
    
    /**
     *  This method is called at the beginning of Solve(). Subclass assemblers can 
     *  use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()
    {}
    
    
    /**
     * Specify what type of basis functions to use.
     * 
     * @param pBasisFunction Basis function to use for normal elements.
     * @param pSurfaceBasisFunction Basis function to use for boundary elements.
     */
    void SetBasisFunctions(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                           AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction)
    {
        if (mWeAllocatedBasisFunctionMemory)
        {
            delete mpBasisFunction;
            delete mpSurfaceBasisFunction;
            mWeAllocatedBasisFunctionMemory = false;
        }
        mpBasisFunction = pBasisFunction;
        mpSurfaceBasisFunction = pSurfaceBasisFunction;
    }
    
    
    /**
     * Set the number of quadrature points to use, per dimension.
     * 
     * This method will throw an exception if the requested number of quadrature
     * points is not supported. (///\todo: There may be a small memory leak if this
     * occurs.)
     * 
     * @param numQuadPoints Number of quadrature points to use per dimension.
     */
    void SetNumberOfQuadraturePointsPerDimension(int numQuadPoints)
    {
        if (mpQuadRule) delete mpQuadRule;
        mpQuadRule = new GaussianQuadratureRule<ELEMENT_DIM>(numQuadPoints);
        if (mpSurfaceQuadRule) delete mpSurfaceQuadRule;
        mpSurfaceQuadRule = new GaussianQuadratureRule<ELEMENT_DIM-1>(numQuadPoints);
    }
    
    
    /**
     * Set the mesh.
     */
    void SetMesh(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    {
        mpMesh = pMesh;
    }
    
    
    /**
     * Set the boundary conditions.
     */
    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions)
    {
        mpBoundaryConditions = pBoundaryConditions;
    }
    
    
};
#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
