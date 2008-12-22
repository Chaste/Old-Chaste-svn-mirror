/*

Copyright (C) University of Oxford, 2008

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


#ifndef TESTPAPILLARYFIBRECALCULATOR_HPP_
#define TESTPAPILLARYFIBRECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "SimpleDataWriter.hpp"
#include "UblasCustomFunctions.hpp"

class TestPapillaryFibreCalculator : public CxxTest::TestSuite
{
public:

    void TestPapillaryFibre(void) throw(Exception)
    {
        std::cout<<"\n Hello, beginning! \n"<<std::flush;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cylinder_14748_elem");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        std::ifstream pap_facefile;
        //coordsfile.open("notforrelease/test/data/OxfordHeart_i_triangles/heartT_renum_i.pts");
        //pap_facefile.open("mesh/test/data/cylinder_muscle.surflist");
        std::cout<<"\n Hello! \n"<<std::flush;
        
        /////////////////////////////////////////////////////////////
        // Defines the numbers of nodes and face nodes
        /////////////////////////////////////////////////////////////
        unsigned int num_elements = mesh.GetNumElements();
   
        // Defines a list into which radial vectors are stored
        c_vector<double, 3> gradients;

        double nearest_r_squared=DBL_MAX,r_squared;
        int nearest_face_node = 0;
        
        // Loops over all elements finding radius vector
        TetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
        while (iter != mesh.GetElementIteratorEnd())
        {
            c_vector<double, 3> centroid = (*iter)->CalculateCentroid();
            // Loops over all papillary face nodes
            c_vector<double,3> coordinates;
            
            
//            for (unsigned int j=0;j<num_pap_face;j++)
            TetrahedralMesh<3,3>::BoundaryNodeIterator bound_node_iter = mesh.GetBoundaryNodeIteratorBegin();
            while (bound_node_iter != mesh.GetBoundaryNodeIteratorEnd())
            {
//                for (unsigned k=0;k<3;k++)
//                {
//                    //coord[k]=coords[pap_face[j]][k];
//                    coordinates[k]=(mesh.GetNode(j)->rGetLocation()[k]);
//                }
                unsigned bound_node_index =  (*bound_node_iter)->GetIndex();           
                coordinates=mesh.GetNode(bound_node_index)->rGetLocation();


                // Calculates the distance between the papillary face node and the centroid
                r_squared =  norm_2(centroid-coordinates);
                // Checks to see if it is the smallest so far - if it is, update the current smallest distance
                if (r_squared < nearest_r_squared)
                {
                    nearest_r_squared = r_squared;
                    nearest_face_node = bound_node_index;
                }   
            }
            // Once we have the papillary face node which is closest, use this to re-define its coordinates
//            for (unsigned k=0;k<3;k++)
//            {
//                coordinates[k]=(mesh.GetNode(nearest_face_node)->rGetLocation()[k]);
//            }
            coordinates = mesh.GetNode(nearest_face_node)->rGetLocation();
            c_vector<double,3> gradients = coordinates-centroid;
        }
        
        std::cout << "After first nested loop" << std::endl;
        
        // Writes-out the radius vector file, we copy the c vector into an std vector to be used by SimpleDataWriter
        std::vector<double> copy_for_writing(3);
        for (unsigned int index=0;index<3;index++)
        {
            copy_for_writing[index]=gradients[index];
        }  
        SimpleDataWriter writer("Fibres", "radius_vector.dat",copy_for_writing, false);
        
        /////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////
        // This part of the code now uses the radius vectors calculated above
        // to generate a structure tensor for each papillary element. The eigen
        // vector with the smallest eigen value corresponds to the axial direction
        // of the muslce
        //////////////////////////////////////////////////////////////////
          
        // Defines entries in (unsmoothed) structure tensor component vectors, where each entry of [I11,I12,I13;I21,I22,I23;I31,I32,I33] is defined to be
        // g^T * g, where g = (g1,g2,g3) is the x,y,z components of the radius vector at that element.
        std::vector< c_matrix<double,3,3> > tensorI(num_elements, zero_matrix<double>(3,3));
        std::vector< c_matrix<double,3,3> > tensorS(num_elements, zero_matrix<double>(3,3));
        
        // Assigns entries in (unsmoothed) structure tensor component vectors, where each entry of [I11,I12,I13;I21,I22,I23;I31,I32,I33] is defined to be
        // g^T * g, where g = (g1,g2,g3) is the x,y,z components of the radius vector at that element.
       // int n,m;     
        for (unsigned int i=0;i<num_elements;i++)
        {       
            c_vector<double, 3> gradient;
            for (unsigned j=0;j<3;j++)
            {
                gradient[j]=gradients[3*i+j];
            }
            tensorI[i] = outer_prod(gradient, gradient);
        }   
        ////////////////////////////////////////////////////////////////////////////////////////////
        // Smoothes the structure tensor components for each papillary elements by looping over all other papillary elements, calculating
        // distance geometric distance between the two elements; if it is within a certain limit, include this in the Gaussian kernel
        ////////////////////////////////////////////////////////////////////////////////////////////
        double g_factor = 0;
        double sigma = 0.5;
        double g_factor_sum = 0;
        double r_max = 1.0;
   
        TetrahedralMesh<3,3>::ElementIterator elem_iter = mesh.GetElementIteratorBegin();
        while (elem_iter != mesh.GetElementIteratorEnd())
        {
            c_vector<double, 3> centroid = (*iter)->CalculateCentroid();  
            g_factor_sum = 0;
            TetrahedralMesh<3,3>::ElementIterator iter_2 = mesh.GetElementIteratorBegin();
            while (iter_2 != mesh.GetElementIteratorEnd())
            {
            c_vector<double, 3> centroid_2 = (*iter_2)->CalculateCentroid();
            double r = norm_2(centroid-centroid_2);             
                if (r < r_max)
                {
                    g_factor = exp(-r/(2*sigma*sigma));
                    
                    g_factor_sum = g_factor + g_factor_sum;
                    
                    for (int l=0;l<3;l++)
                    {
                        for (int k=0;k<3;k++)
                        {
                            tensorS[ (*iter)->GetIndex()](k,l) = tensorS[ (*iter)->GetIndex()](k,l) + g_factor*tensorI[ (*iter_2)->GetIndex()](k,l);
                        }
                    }     
                }
        }         
        for (unsigned int l=0;l<3;l++)
        {
            for (unsigned int k=0;k<3;k++)
            {
                tensorS[ (*iter)->GetIndex()](k,l) = tensorS[ (*iter)->GetIndex()](k,l)/g_factor_sum;
            }
        } 
        //calculate Eigenvector for smalles eigen value
        c_vector<double, 3> eigenvector;
        c_matrix<double, 3, 3> A;
        A= tensorS[ (*iter)->GetIndex()];   
        eigenvector = CalculateEigenvectorForSmallestEigenvalue(A);        
        }      
    }/*end of the test*/
};

#endif /*TESTPAPILLARYFIBRECALCULATOR_HPP_*/
