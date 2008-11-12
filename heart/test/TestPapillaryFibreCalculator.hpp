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

std::string epi_face_file = "/home/chaste/heart_data/pap_face_n.tri";

ifstream coordsfile;
ifstream pap_facefile;
coordsfile.open("/home/chaste/heart_data/heartT_renum_i.pts");
pap_facefile.open("/home/chaste/heart_data/pap_face_n.dat");


/////////////////////////////////////////////////////////////
  // Defines the numbers of nodes and face nodes
  /////////////////////////////////////////////////////////////
  int num_nodes = 4310704;
  int num_elements = 24217344;
    int num_pap_face = 42539;
    


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
  ifstream centsfile("/home/chaste/heart_data/centroids.pts");
 
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
ifstream papselemsfile("/home/chaste/heart_data/paps_elems.dat");
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
 
 

  // Reads in coords file
  double x,y,z,dummy;
  coordsfile >> dummy;
  for(int i=0;i<num_nodes;i++)
    {
      coordsfile >> x >> y >> z;
      coords[i][0] = x*x_factor;
      coords[i][1] = y*y_factor;
      coords[i][2] = z*z_factor;

    }
cout << coords[10][0] << " " << coords[10][2] << "\n ";


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
		
		
		 // Writes-out the radius vector file
  ofstream vectorfile("/home/chaste/heart_data/radius_vector.dat");
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
 double *I11;
  I11 = new double [num_paps_elems];

double *I12;
  I12 = new double [num_paps_elems];

double *I13;
  I13 = new double [num_paps_elems];

double *I21;
  I21 = new double [num_paps_elems];

double *I22;
  I22 = new double [num_paps_elems];

double *I23;
  I23 = new double [num_paps_elems];

double *I31;
  I31 = new double [num_paps_elems];

double *I32;
  I32 = new double [num_paps_elems];

double *I33;
  I33 = new double [num_paps_elems];
  
  double *I11s;
  I11s = new double [num_paps_elems];

double *I12s;
  I12s = new double [num_paps_elems];

double *I13s;
  I13s = new double [num_paps_elems];

double *I21s;
  I21s = new double [num_paps_elems];

double *I22s;
  I22s = new double [num_paps_elems];

double *I23s;
  I23s = new double [num_paps_elems];

double *I31s;
  I31s = new double [num_paps_elems];

double *I32s;
  I32s = new double [num_paps_elems];

double *I33s;
  I33s = new double [num_paps_elems];

// Initialises structure tensor component vectors
  for(int i=0;i<num_paps_elems;i++)
    {
      I11[i] = 0;
      I12[i] = 0;
      I13[i] = 0;
      I21[i] = 0;
      I22[i] = 0;
      I23[i] = 0;
      I31[i] = 0;
      I32[i] = 0;
      I33[i] = 0;
      
      I11s[i] = 0;
      I12s[i] = 0;
      I13s[i] = 0;
      I21s[i] = 0;
      I22s[i] = 0;
      I23s[i] = 0;
      I31s[i] = 0;
      I32s[i] = 0;
      I33s[i] = 0;
    }
 
// Assigns entries in (unsmoothed) structure tensor component vectors, where each entry of [I11,I12,I13;I21,I22,I23;I31,I32,I33] is defined to be
// g^T * g, where g = (g1,g2,g3) is the x,y,z components of the radius vector at that element. 
 int n,m;

  for(int i=0;i<num_paps_elems;i++)
    {
      n = paps_elems[i];
      
      c_vector<double, 3> gradient;
      for (unsigned j=0;j<3;j++)
      {
        gradient[j]=gradients[n][j];
      }
      
      c_matrix <double, 3, 3> tensor = outer_prod(gradient, gradient);
      
      I11[i] = tensor(0,0);
      I12[i] = tensor(0,1);
      I13[i] = tensor(0,2);
      I21[i] = tensor(1,0);
      I22[i] = tensor(1,1);
      I23[i] = tensor(1,2);
      I31[i] = tensor(2,0);
      I32[i] = tensor(2,1);
      I33[i] = tensor(2,2);
      
      TS_ASSERT_DELTA(I11[i], tensor(0,0), 1e-16);    

    }
 
 ////////////////////////////////////////////////////////////////////////////////////////////
 // Smoothes the structure tensor components for each papillary elements by looping over all other papillary elements, calculating
 // distance geometric distance between the two elements; if it is within a certain limit, include this in the Gaussian kernel
 ////////////////////////////////////////////////////////////////////////////////////////////
  double x_p,y_p,z_p,sv_I11,sv_I12,sv_I13,sv_I21,sv_I22,sv_I23,sv_I31,sv_I32,sv_I33,r_max;
  double g_factor = 0;
  double sigma = 0.5;
  double g_factor_sum = 0;
  r_max = 1.0;

  for(int i=0;i<num_paps_elems;i++)
    {

      if(i == 100 || i == 1000 || i == 10000 || i == 100000 || i == 500000 || i == 800000 || i == 1000000)
        cout << i << "\n";
        
      sv_I11 = 0;
      sv_I12 = 0;
      sv_I13 = 0;
      sv_I21 = 0;
      sv_I22 = 0;
      sv_I23 = 0;
      sv_I31 = 0;
      sv_I32 = 0;
      sv_I33 = 0;

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
    
              sv_I11 = sv_I11 + g_factor*I11[j];
              sv_I12 = sv_I12 + g_factor*I12[j];
              sv_I13 = sv_I13 + g_factor*I13[j];
              sv_I21 = sv_I21 + g_factor*I21[j];
              sv_I22 = sv_I22 + g_factor*I22[j];
              sv_I23 = sv_I23 + g_factor*I23[j];
              sv_I31 = sv_I31 + g_factor*I31[j];
              sv_I32 = sv_I32 + g_factor*I32[j];
              sv_I33 = sv_I33 + g_factor*I33[j];

            }
        }

      I11s[i] = sv_I11/g_factor_sum;
      I12s[i] = sv_I12/g_factor_sum;
      I13s[i] = sv_I13/g_factor_sum;
      I21s[i] = sv_I21/g_factor_sum;
      I22s[i] = sv_I22/g_factor_sum;
      I23s[i] = sv_I23/g_factor_sum;
      I31s[i] = sv_I31/g_factor_sum;
      I32s[i] = sv_I32/g_factor_sum;
      I33s[i] = sv_I33/g_factor_sum;

    }
    
    // Writes-out the list of nearest neighbours
  ofstream I11sfile("/home/chaste/heart_data/I11s.dat");
  ofstream I12sfile("/home/chaste/heart_data/I12s.dat");
  ofstream I13sfile("/home/chaste/heart_data/I13s.dat");
  ofstream I21sfile("/home/chaste/heart_data/I21s.dat");
  ofstream I22sfile("/home/chaste/heart_data/I22s.dat");
  ofstream I23sfile("/home/chaste/heart_data/I23s.dat");
  ofstream I31sfile("/home/chaste/heart_data/I31s.dat");
  ofstream I32sfile("/home/chaste/heart_data/I32s.dat");
  ofstream I33sfile("/home/chaste/heart_data/I33s.dat");
  
  for(int i=0;i<num_paps_elems;i++)
    {
     I11sfile << I11s[i] << "\n";
     I12sfile << I12s[i] << "\n";
     I13sfile << I13s[i] << "\n";
     I21sfile << I21s[i] << "\n";
     I22sfile << I22s[i] << "\n";
     I23sfile << I23s[i] << "\n";
     I31sfile << I31s[i] << "\n";
     I32sfile << I32s[i] << "\n";
     I33sfile << I33s[i] << "\n";

    }
    
	    coordsfile.close();
	    pap_facefile.close();
        

    }
};

#endif /*TESTPAPILLARYFIBRECALCULATOR_HPP_*/
