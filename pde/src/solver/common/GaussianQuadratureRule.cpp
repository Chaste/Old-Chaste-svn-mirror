/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include <cmath>

#include "GaussianQuadratureRule.hpp"
#include "Exception.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned ELEMENT_DIM>
const ChastePoint<ELEMENT_DIM>& GaussianQuadratureRule<ELEMENT_DIM>::rGetQuadPoint(unsigned index) const
{
    assert(index < mNumQuadPoints);
    return mPoints[index];
}

template<unsigned ELEMENT_DIM>
double GaussianQuadratureRule<ELEMENT_DIM>::GetWeight(unsigned index) const
{
    assert(index < mNumQuadPoints);
    return mWeights[index];
}

template<unsigned ELEMENT_DIM>
unsigned GaussianQuadratureRule<ELEMENT_DIM>::GetNumQuadPoints() const
{
    return mNumQuadPoints;
}

template<unsigned ELEMENT_DIM>
GaussianQuadratureRule<ELEMENT_DIM>::GaussianQuadratureRule(unsigned numPointsInEachDimension)
{
    mNumQuadPoints = (unsigned) SmallPow((double) numPointsInEachDimension,(ELEMENT_DIM));

    mWeights.reserve(mNumQuadPoints);
    mPoints.reserve(mNumQuadPoints);

    switch (ELEMENT_DIM)
    {
        case 0 :
        {
            // mNumQuadPoints should have been set to be one as
            // it is numPointsInEachDim^{0}
            mWeights.push_back(1);
            mPoints.push_back(ChastePoint<ELEMENT_DIM>());
        }
        break;
        case 1 :
        {
            switch (numPointsInEachDimension)
            {
                case 1: // 1d, 1 point
                    mWeights.push_back(1);
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.5)); //check
                    break;

                case 2: // 1d, 2 points
                    mWeights.push_back(0.5);
                    mWeights.push_back(0.5);

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.21132486540519));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.78867513459481));
                    break;

                case 3: // 1d, 3 points
                    mWeights.push_back(5.0/18.0);
                    mWeights.push_back(4.0/9.0);
                    mWeights.push_back(5.0/18.0);

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.1127016654));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.5));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.8872983346));
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
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.25,0.5));
                    break;

                case 2: // 2d, 2 points per dimension
                    mWeights.push_back(0.19716878364870);
                    mWeights.push_back(0.19716878364870);
                    mWeights.push_back(0.05283121635130);
                    mWeights.push_back(0.05283121635130);

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.16666666666667,0.21132486540519));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62200846792815,0.21132486540519));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04465819873852,0.78867513459481));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.16666666666667,0.78867513459481));
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

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,0.11270166540000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,0.11270166540000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.78729833458393,0.11270166540000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,0.50000000000000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.25000000000000,0.50000000000000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,0.50000000000000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.01270166538393,0.88729833460000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,0.88729833460000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,0.88729833460000));
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
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.25000000000000,0.50000000000000,0.12500000000000));
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

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.16666666666667,   0.21132486540519,   0.13144585576580));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62200846792815,   0.21132486540519,   0.03522081090086));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04465819873852,   0.78867513459481,   0.03522081090086));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.16666666666667,   0.78867513459481,   0.00943738783766));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.16666666666667,   0.21132486540519,   0.49056261216234));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62200846792815,   0.21132486540519,   0.13144585576580));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04465819873852,   0.78867513459481,   0.13144585576580));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.16666666666667,   0.78867513459481,   0.03522081090086));
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

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,   0.11270166540000,   0.08872983347426));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,   0.11270166540000,   0.05000000000803));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.78729833458393,   0.11270166540000,   0.01127016654181));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,   0.50000000000000,   0.05000000000803));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.25000000000000,   0.50000000000000,   0.02817541635000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,   0.50000000000000,   0.00635083269197));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.01270166538393,   0.88729833460000,   0.01127016654181));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,   0.88729833460000,   0.00635083269197));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,   0.88729833460000,   0.00143149884212));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,   0.11270166540000,   0.39364916729197));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,   0.11270166540000,   0.22182458365000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.78729833458393,   0.11270166540000,   0.05000000000803));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,   0.50000000000000,   0.22182458365000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.25000000000000,   0.50000000000000,   0.12500000000000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,   0.50000000000000,   0.02817541635000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.01270166538393,   0.88729833460000,   0.05000000000803));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,   0.88729833460000,   0.02817541635000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,   0.88729833460000,   0.00635083269197));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,   0.11270166540000,   0.69856850110968));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,   0.11270166540000,   0.39364916729197));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.78729833458393,   0.11270166540000,   0.08872983347426));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,   0.50000000000000,   0.39364916729197));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.25000000000000,   0.50000000000000,   0.22182458365000));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44364916730000,   0.50000000000000,   0.05000000000803));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.01270166538393,   0.88729833460000,   0.08872983347426));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.05635083270000,   0.88729833460000,   0.05000000000803));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10000000001607,   0.88729833460000,   0.01127016654181));
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

                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.06943184420000,   0.06012499793653));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.06943184420000,   0.04328879995478));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.06943184420000,   0.02132226325621));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.86595709258901,   0.06943184420000,   0.00448606527446));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.33000947820000,   0.04328879995478));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.33000947820000,   0.03116707302848));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44888729930184,   0.33000947820000,   0.01535160449661));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.33000947820000,   0.00322987757031));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.66999052180000,   0.02132226325621));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10890625570184,   0.66999052180000,   0.01535160449661));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.66999052180000,   0.00756156217830));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.66999052180000,   0.00159090341870));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.00482078098901,   0.93056815580000,   0.00448606527446));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.93056815580000,   0.00322987757031));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.93056815580000,   0.00159090341870));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.93056815580000,   0.00033471571455));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.06943184420000,   0.28577404826889));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.06943184420000,   0.20575161800155));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.06943184420000,   0.10134469352354));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.86595709258901,   0.06943184420000,   0.02132226325621));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.33000947820000,   0.20575161800155));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.33000947820000,   0.14813706341321));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44888729930184,   0.33000947820000,   0.07296615908496));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.33000947820000,   0.01535160449661));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.66999052180000,   0.10134469352354));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10890625570184,   0.66999052180000,   0.07296615908496));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.66999052180000,   0.03594009661688));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.66999052180000,   0.00756156217830));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.00482078098901,   0.93056815580000,   0.02132226325621));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.93056815580000,   0.01535160449661));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.93056815580000,   0.00756156217830));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.93056815580000,   0.00159090341870));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.06943184420000,   0.58018304432012));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.06943184420000,   0.41772022627335));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.06943184420000,   0.20575161800155));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.86595709258901,   0.06943184420000,   0.04328879995478));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.33000947820000,   0.41772022627335));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.33000947820000,   0.30075023588863));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44888729930184,   0.33000947820000,   0.14813706341321));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.33000947820000,   0.03116707302848));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.66999052180000,   0.20575161800155));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10890625570184,   0.66999052180000,   0.14813706341321));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.66999052180000,   0.07296615908496));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.66999052180000,   0.01535160449661));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.00482078098901,   0.93056815580000,   0.04328879995478));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.93056815580000,   0.03116707302848));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.93056815580000,   0.01535160449661));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.93056815580000,   0.00322987757031));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.06943184420000,   0.80583209465249));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.06943184420000,   0.58018304432012));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.06943184420000,   0.28577404826889));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.86595709258901,   0.06943184420000,   0.06012499793653));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.33000947820000,   0.58018304432012));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.33000947820000,   0.41772022627335));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.44888729930184,   0.33000947820000,   0.20575161800155));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.62347184427491,   0.33000947820000,   0.04328879995478));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.66999052180000,   0.28577404826889));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.10890625570184,   0.66999052180000,   0.20575161800155));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.22110322249816,   0.66999052180000,   0.10134469352354));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.30709631152509,   0.66999052180000,   0.02132226325621));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.00482078098901,   0.93056815580000,   0.06012499793653));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.02291316667491,   0.93056815580000,   0.04328879995478));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.04651867752509,   0.93056815580000,   0.02132226325621));
                    mPoints.push_back(ChastePoint<ELEMENT_DIM>(0.06461106321099,   0.93056815580000,   0.00448606527446));
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


//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class GaussianQuadratureRule<0>;
template class GaussianQuadratureRule<1>;
template class GaussianQuadratureRule<2>;
template class GaussianQuadratureRule<3>;
template class GaussianQuadratureRule<4>;
