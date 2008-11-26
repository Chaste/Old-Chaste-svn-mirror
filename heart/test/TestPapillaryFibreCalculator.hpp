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
using namespace std;
class TestPapillaryFibreCalculator : public CxxTest::TestSuite
{
public:

    void TestPapillaryFibre(void) throw(Exception)
    {

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cylinder_muscle");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

//ifstream coordsfile;
ifstream pap_facefile;
//coordsfile.open("notforrelease/test/data/OxfordHeart_i_triangles/heartT_renum_i.pts");
pap_facefile.open("mesh/test/data/cylinder_muscle.surflist");


/////////////////////////////////////////////////////////////
  // Defines the numbers of nodes and face nodes
  /////////////////////////////////////////////////////////////
  int num_nodes = mesh.GetNumNodes();
  int num_elements = mesh.GetNumElements();
  int num_pap_face = mesh.GetNumBoundaryNodes();
    


 // Defines coordinate list 
  double **coords;
  coords = new double *[num_nodes];
  for(int i=0;i<num_nodes;i++)
    coords[i] = new double[3];

// Defines list of all papillary face nodes
  int *pap_face;
  pap_face = new int [num_pap_face];
/////////////////////////////////////////////////////////////////
  // Reads the centroids of the elements
  /////////////////////////////////////////////////////////////////
  ifstream centsfile("notforrelease/test/data/OxfordHeart_i_triangles/centroids.pts");
 
  // Defines the number of centroids (same as number of elements...)
  int num_cents = 24217344;
 // Defines scaling factors for coords file
  double x_factor = 0.053/1000;
  double y_factor = x_factor;
  double z_factor = 0.049/1000;
  
  double **cents;
  cents = new double *[num_cents];
  for(int i=0;i<num_cents;i++)
    cents[i] = new double[3];
  
  std::cout << "b"<<"\n";
 
  // Reads in cents file
    for(int i=0;i<num_cents;i++)
    {
      double x,y,z;
      centsfile >> x >> y >> z;
      cents[i][0] = x*x_factor;
      cents[i][1] = y*y_factor;
      cents[i][2] = z*z_factor;
    }
cout << cents[10][0] << " " << cents[10][2] << "\n ";
  
// Reads-in list of only papillary elements
ifstream papselemsfile("notforrelease/test/data/OxfordHeart_i_triangles/paps_elems.dat");
int num_paps_elems = 1278265;
int *paps_elems;
  paps_elems = new int [num_paps_elems];

int pap_element;
for(int i=0;i<num_paps_elems;i++)
    {
      papselemsfile >> pap_element;
      paps_elems[i] = pap_element;
    }
 cout << paps_elems[1278262] << "\n";
 
 

   // Reads in list of papillary nodes
  int pap_face_value;
for(int i=0;i<num_pap_face;i++)
    {
      pap_facefile >> pap_face_value;
      pap_face[i] = pap_face_value;
    }
cout << pap_face[10] << "\n ";

// Defines a list into which radial vectors are stored
double **gradients;
  gradients = new double *[num_elements];
  for(int i=0;i<num_elements;i++)
    gradients[i] = new double[3];

  // Defines quantities used below
  double x_c,y_c,z_c,x_f,y_f,z_f;
  double nearest_r_squared,r_squared;
  int nearest_face_node = 0;

		
		 // Loops over all elements finding radius vector
//  for(int i=0;i<num_paps_elems;i++)
    for(int i=0;i<10;i++)
    {
            // Sets the distance to be very big by default
          nearest_r_squared = DBL_MAX;
          c_vector<double,3> centroid;
          for (unsigned j=0;j<3;j++)
          {
            centroid[j]=cents[paps_elems[i]][j];
          }
          // Loops over all papillary face nodes
          for(int j=0;j<num_pap_face;j++)
            {
              // Defines the coordinates of the papillary face node
              x_f = coords[pap_face[j]][0];
              y_f = coords[pap_face[j]][1];
              z_f = coords[pap_face[j]][2];
              c_vector<double,3> coord;
              for (unsigned k=0;k<3;k++)
              {
                coord[k]=coords[pap_face[j]][k];
              }
     
              // Calculates the distance between the papillary face node and the centroid
              r_squared =  norm_2(centroid-coord);

              // Checks to see if it is the smallest so far - if it is, update the current smallest distance
              if(r_squared < nearest_r_squared)
                {
                  nearest_r_squared = r_squared;
                  nearest_face_node = j;
                }

            }
          // Once we have the papillary face node which is closest, use this to re-define its coordinates
          x_f = coords[pap_face[nearest_face_node]][0];
          y_f = coords[pap_face[nearest_face_node]][1];
          z_f = coords[pap_face[nearest_face_node]][2];

          // Defines the radius vector to be r = r1 - r2 
          gradients[i][0] = x_c - x_f;
          gradients[i][1] = y_c - y_f;
          gradients[i][2] = z_c - z_f;

 
    }
		
std::cout << "c"<<"\n";        
		
		 // Writes-out the radius vector file
  ofstream vectorfile("notforrelease/test/data/OxfordHeart_i_triangles/radius_vector.dat");
  for(int i=0;i<num_elements;i++)
    {
      vectorfile << gradients[i][0] << " " << gradients[i][1] << " " << gradients[i][2] << "\n";

    }
  vectorfile.close();

  cout << "done! \n";
	
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    // This part of the code now uses the radius vectors calculated above
    // to generate a structure tensor for each papillary element. The eigen
    // vector with the smallest eigen value corresponds to the axial direction
    // of the muslce
    //////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    
    
  
 // Defines entries in (unsmoothed) structure tensor component vectors, where each entry of [I11,I12,I13;I21,I22,I23;I31,I32,I33] is defined to be
// g^T * g, where g = (g1,g2,g3) is the x,y,z components of the radius vector at that element. 

     std::vector< c_matrix<double,3,3> > tensorI(num_paps_elems, zero_matrix<double>(3,3));
     std::vector< c_matrix<double,3,3> > tensorS(num_paps_elems, zero_matrix<double>(3,3));     
//     c_matrix<double,3,3> computed(zero_matrix<double>(3,3));
     
// Initialises structure tensor component vectors

//     for(int k=0;k<3;k++)
//     {
//        for(int j=0;j<3;j++)
//        {
//         computed(j,k)=0;
//        }  
//     }
     
//     for(int i=0;i<num_paps_elems;i++)
//     {
//        tensorI[i]=computed;
//        tensorS[i]=computed;
//     }
    //TS_ASSERT_EQUALS(tensorI[1](0,0), 24217344U);
    
  

// Assigns entries in (unsmoothed) structure tensor component vectors, where each entry of [I11,I12,I13;I21,I22,I23;I31,I32,I33] is defined to be
// g^T * g, where g = (g1,g2,g3) is the x,y,z components of the radius vector at that element. 
 int n,m;

  for(int i=0;i<10;i++)
    {
      n = paps_elems[i];
      
      c_vector<double, 3> gradient;
      for (unsigned j=0;j<3;j++)
      {
        gradient[j]=gradients[n][j];
      }
      
//      c_matrix <double, 3, 3> tensor = outer_prod(gradient, gradient);
      
//      for(int k=0;k<3;k++)
//      {
//        for(int j=0;j<3;j++)
//        {
//            tensorI[i](j,k) = tensor(j,k);
//        }
//      }
      tensorI[i] = outer_prod(gradient, gradient);  
      
      
//      TS_ASSERT_DELTA(I11[i], tensor(0,0), 1e-16);    

    }
 
//   assert(0);
std::cout << "d"<<"\n";
 
 ////////////////////////////////////////////////////////////////////////////////////////////
 // Smoothes the structure tensor components for each papillary elements by looping over all other papillary elements, calculating
 // distance geometric distance between the two elements; if it is within a certain limit, include this in the Gaussian kernel
 ////////////////////////////////////////////////////////////////////////////////////////////
  double x_p,y_p,z_p,r_max;
  double g_factor = 0;
  double sigma = 0.5;
  double g_factor_sum = 0;
  r_max = 1.0;

      

  for(int i=0;i<10;i++)
    {

      if(i == 100 || i == 1000 || i == 10000 || i == 100000 || i == 500000 || i == 800000 || i == 1000000)
        cout << i << "\n";
        
      g_factor_sum = 0;

      n = paps_elems[i];

      x_p = cents[n][0];
      y_p = cents[n][1];
      z_p = cents[n][2];

      for(int j=0;j<num_paps_elems;j++)
        {
          m = paps_elems[j];

          x_c = cents[m][0];
          y_c = cents[m][1];
          z_c = cents[m][2];

          double r = sqrt((x_p - x_c)*(x_p - x_c) + (y_p - y_c)*(y_p - y_c) + (z_p - z_c)*(z_p - z_c) );

          if(r < r_max)
            {
              g_factor = exp(-( (x_c - x_p)*(x_c - x_p) + (y_c - y_p)*(y_c - y_p) + (z_c - z_p)*(z_c - z_p) )/(2*sigma*sigma));

              g_factor_sum = g_factor + g_factor_sum;
    
              for(int l=0;l<3;l++)
              {
                for(int k=0;k<3;k++)
                {
                    tensorS[i](k,l) = tensorS[i](k,l) + g_factor*tensorI[j](k,l);
                }
              }
             

            }
        }
       
       for(int l=0;l<3;l++)
       {
           for(int k=0;k<3;k++)
           {
               tensorS[i](k,l) = tensorS[i](k,l)/g_factor_sum;
           }
       }

    }
    
    
	    pap_facefile.close();
        

    }
};

#endif /*TESTPAPILLARYFIBRECALCULATOR_HPP_*/
