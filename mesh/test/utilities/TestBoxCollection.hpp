/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef TESTBOXCOLLECTION_HPP_
#define TESTBOXCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "TetrahedralMesh.hpp"
#include "BoxCollection.hpp"
#include "TrianglesMeshReader.hpp"
#include "Debug.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBoxCollection : public CxxTest::TestSuite
{
public:
    void TestBox() throw (Exception)
    {
			c_vector<double, 2*2> box_size;
			box_size(0) = -0.1; // min x
			box_size(1) = 1.1; // max x
			box_size(2) = -0.1; // min y
			box_size(3) = 1.1; // max y

			Box<2> test_box(box_size);
			c_vector<double, 2*2> returned_min_max_values = test_box.rGetMinAndMaxValues();
			for (unsigned i=0; i<4; i++)
			{
				TS_ASSERT_EQUALS(returned_min_max_values(i), box_size(i));
			}

			c_vector<double, 2> node_location;
			node_location(0) = 0.5;
			node_location(1) = 0.5;

			Node<2> test_node(213, node_location);

			test_box.AddNode(&test_node);
			std::set< Node<2>* > nodes_contained_before = test_box.rGetNodesContained();

			TS_ASSERT_EQUALS(*(nodes_contained_before.begin()), &test_node);
			TS_ASSERT_EQUALS((*(nodes_contained_before.begin()))->GetIndex(), 213u);

			test_box.RemoveNode(&test_node);
			std::set< Node<2>* > nodes_contained_after = test_box.rGetNodesContained();
			TS_ASSERT(nodes_contained_after.empty());

    }


    void TestBoxGeneration1d() throw (Exception)
    {
    		// Create a mesh
			TetrahedralMesh<1,1> mesh;
			mesh.ConstructLinearMesh(20);

			double cut_off_length = 5.0;
			c_vector<double, 2> domain_size;
			domain_size(0) = -0.1;
			domain_size(1) = 20.15;

			BoxCollection<1> box_collection(cut_off_length, domain_size);

			if(PetscTools::IsSequential())
			{
				box_collection.SetupAllLocalBoxes();
				TS_ASSERT_THROWS_THIS(box_collection.SetupAllLocalBoxes(), "Local Boxes Already Set");

				for (unsigned i=0; i<mesh.GetNumNodes(); i++)
				{
					unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
					// If we own the box, add the node.
					if(box_collection.GetBoxOwnership(box_index))
					{
						box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
					}
				}

				TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 5u);

				for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
				{
					if(box_collection.GetBoxOwnership(i))
					{
						std::set< Node<1>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
						c_vector<double, 2> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

						for (std::set< Node<1>* >::iterator it_nodes_in_box = nodes_in_box.begin();
							 it_nodes_in_box != nodes_in_box.end();
							 it_nodes_in_box++)
						{
							Node<1>* current_node = *it_nodes_in_box;
							double x_position = current_node->rGetLocation()[0];

							double epsilon = 1e-12;

							TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
							TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
						}
					}
				}

				if(box_collection.GetBoxOwnership(0))
				{
					std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
					std::set<unsigned> correct_answer_0;
					correct_answer_0.insert(0);
					correct_answer_0.insert(1);
					TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
				}
				if(box_collection.GetBoxOwnership(1))
				{
					std::set<unsigned> local_boxes_to_box_1 = box_collection.GetLocalBoxes(1);
					std::set<unsigned> correct_answer_1;
					correct_answer_1.insert(0);
					correct_answer_1.insert(1);
					correct_answer_1.insert(2);
					TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);
				}
				if(box_collection.GetBoxOwnership(4))
				{
					std::set<unsigned> local_boxes_to_box_4 = box_collection.GetLocalBoxes(4);
					std::set<unsigned> correct_answer_4;
					correct_answer_4.insert(3);
					correct_answer_4.insert(4);
					TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);
				}

				c_vector<double,1> miles_away;
				miles_away(0) = 47323854;
				TS_ASSERT_THROWS_CONTAINS(box_collection.CalculateContainingBox(miles_away), "The point provided is outside all of the boxes");
			}
			else
			{
				TS_ASSERT_THROWS_CONTAINS(box_collection.SetupAllLocalBoxes(), "Using SetupAllLocalBoxes() in parallel. This will likely lead to errors. If possible use SetupLocalBoxesHalfOnly()");

			}

			// Coverage
			c_vector<unsigned, 1> indices=box_collection.CalculateCoordinateIndices(0);
			TS_ASSERT_EQUALS(indices[0],0u);
    }

    // very simple test
    void TestAddElement() throw(Exception)
    {
		TetrahedralMesh<1,1> mesh;
		mesh.ConstructRegularSlabMesh(0.5, 1.0);

		double width = 0.4;
		c_vector<double, 2> domain_size;
		domain_size(0) = 0.0;
		domain_size(1) = 1.0;

		BoxCollection<1> box_collection(width, domain_size);
		if(box_collection.GetBoxOwnership(0))
		{
			box_collection.rGetBox(0).AddElement(mesh.GetElement(0));
			TS_ASSERT_EQUALS(box_collection.rGetBox(0).rGetElementsContained().size(), 1u);
		}
		if(box_collection.GetBoxOwnership(1))
		{
			TS_ASSERT_EQUALS(box_collection.rGetBox(1).rGetElementsContained().size(), 0u);
		}
		if(box_collection.GetBoxOwnership(2))
		{
		TS_ASSERT_EQUALS(box_collection.rGetBox(2).rGetElementsContained().size(), 0u);
		}
		if(box_collection.GetBoxOwnership(0))
		{
		TS_ASSERT_EQUALS(*(box_collection.rGetBox(0).rGetElementsContained().begin()), mesh.GetElement(0));
		}

    }


    void TestSetupAllLocalBoxes2d() throw(Exception)
    {
    	// Only in serial
    	if(PetscTools::IsSequential())
    	{
			double width = 1.0;

			c_vector<double, 2*2> domain_size;
			domain_size(0) = 0;
			domain_size(1) = 4-0.01;
			domain_size(2) = 0;
			domain_size(3) = 3-0.01;

			BoxCollection<2> box_collection(width, domain_size);

			assert(box_collection.GetNumBoxes()==12); // 4 * 3 boxes altogether

			box_collection.SetupAllLocalBoxes();

			if(box_collection.GetBoxOwnership(0))
			{
				std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

				std::set<unsigned> correct_answer_0;
				correct_answer_0.insert(0);
				correct_answer_0.insert(1);
				correct_answer_0.insert(4);
				correct_answer_0.insert(5);
				TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
			}

			if(box_collection.GetBoxOwnership(3))
			{
				std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
				std::set<unsigned> correct_answer_3;
				correct_answer_3.insert(3);
				correct_answer_3.insert(2);
				correct_answer_3.insert(6);
				correct_answer_3.insert(7);
				//TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);
			}

			if(box_collection.GetBoxOwnership(5))
			{
				std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
				std::set<unsigned> correct_answer_5;
				correct_answer_5.insert(0);
				correct_answer_5.insert(1);
				correct_answer_5.insert(2);
				correct_answer_5.insert(4);
				correct_answer_5.insert(5);
				correct_answer_5.insert(6);
				correct_answer_5.insert(8);
				correct_answer_5.insert(9);
				correct_answer_5.insert(10);
				TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);
			}

			if(box_collection.GetBoxOwnership(10))
			{
				std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
				std::set<unsigned> correct_answer_10;
				correct_answer_10.insert(5);
				correct_answer_10.insert(6);
				correct_answer_10.insert(7);
				correct_answer_10.insert(9);
				correct_answer_10.insert(10);
				correct_answer_10.insert(11);
				TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);
			}
    	}

    }



    void TestSetupAllLocalBoxes3d() throw(Exception)
    {
    	// Only sequential
    	if(PetscTools::IsSequential())
    	{
			double width = 1.0;

			c_vector<double, 2*3> domain_size;
			domain_size(0) = 0;
			domain_size(1) = 4-0.01;
			domain_size(2) = 0;
			domain_size(3) = 3-0.01;
			domain_size(4) = 0;
			domain_size(5) = 2-0.01;

			BoxCollection<3> box_collection(width, domain_size);

			assert(box_collection.GetNumBoxes()==24); // 4 * 3 * 2 boxes altogether

			box_collection.SetupAllLocalBoxes();

			if(box_collection.GetBoxOwnership(0))
			{
				std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

				std::set<unsigned> correct_answer_0;
				correct_answer_0.insert(0);
				correct_answer_0.insert(1);
				correct_answer_0.insert(4);
				correct_answer_0.insert(5);
				correct_answer_0.insert(12);
				correct_answer_0.insert(13);
				correct_answer_0.insert(16);
				correct_answer_0.insert(17);
				TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
			}

			if(box_collection.GetBoxOwnership(3))
			{
				std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
				std::set<unsigned> correct_answer_3;
				correct_answer_3.insert(3);
				correct_answer_3.insert(2);
				correct_answer_3.insert(6);
				correct_answer_3.insert(7);
				correct_answer_3.insert(14);
				correct_answer_3.insert(15);
				correct_answer_3.insert(18);
				correct_answer_3.insert(19);
				TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);
			}

			if(box_collection.GetBoxOwnership(5))
			{
				std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
				std::set<unsigned> correct_answer_5;
				correct_answer_5.insert(0);
				correct_answer_5.insert(1);
				correct_answer_5.insert(2);
				correct_answer_5.insert(4);
				correct_answer_5.insert(5);
				correct_answer_5.insert(6);
				correct_answer_5.insert(8);
				correct_answer_5.insert(9);
				correct_answer_5.insert(10);
				correct_answer_5.insert(12);
				correct_answer_5.insert(13);
				correct_answer_5.insert(14);
				correct_answer_5.insert(16);
				correct_answer_5.insert(17);
				correct_answer_5.insert(18);
				correct_answer_5.insert(20);
				correct_answer_5.insert(21);
				correct_answer_5.insert(22);

				TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);
			}

			if(box_collection.GetBoxOwnership(19))
			{
				std::set<unsigned> local_boxes_to_box_19 = box_collection.GetLocalBoxes(19);

				std::set<unsigned> correct_answer_19;
				correct_answer_19.insert(2);
				correct_answer_19.insert(3);
				correct_answer_19.insert(6);
				correct_answer_19.insert(7);
				correct_answer_19.insert(10);
				correct_answer_19.insert(11);
				correct_answer_19.insert(14);
				correct_answer_19.insert(15);
				correct_answer_19.insert(18);
				correct_answer_19.insert(19);
				correct_answer_19.insert(22);
				correct_answer_19.insert(23);
				TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);
			}

			if(box_collection.GetBoxOwnership(22))
			{
				std::set<unsigned> local_boxes_to_box_22 = box_collection.GetLocalBoxes(22);
				std::set<unsigned> correct_answer_22;
				correct_answer_22.insert(5);
				correct_answer_22.insert(6);
				correct_answer_22.insert(7);
				correct_answer_22.insert(9);
				correct_answer_22.insert(10);
				correct_answer_22.insert(11);
				correct_answer_22.insert(17);
				correct_answer_22.insert(18);
				correct_answer_22.insert(19);
				correct_answer_22.insert(21);
				correct_answer_22.insert(22);
				correct_answer_22.insert(23);
				TS_ASSERT_EQUALS(local_boxes_to_box_22, correct_answer_22);
			}
    	}
    }

    // Currently only tested in 3d
    // \todo test in 1 and 2d
    void TestSetupLocalBoxesHalfOnly3d() throw(Exception)
	{
		double width = 1.0;

		c_vector<double, 2*3> domain_size;
		domain_size(0) = 0.0;
		domain_size(1) = 3.0;
		domain_size(2) = 0.0;
		domain_size(3) = 3.0;
		domain_size(4) = 0.0;
		domain_size(5) = 2.0;

		BoxCollection<3> box_collection(width, domain_size);

		box_collection.SetupLocalBoxesHalfOnly();

		if(box_collection.GetBoxOwnership(0))
		{
			std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

			std::set<unsigned> correct_answer_0;
			correct_answer_0.insert(0);
			correct_answer_0.insert(1);
			correct_answer_0.insert(3);
			correct_answer_0.insert(4);
			correct_answer_0.insert(9);
			correct_answer_0.insert(10);
			correct_answer_0.insert(12);
			correct_answer_0.insert(13);
			TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
		}
		if(box_collection.GetBoxOwnership(1) && PetscTools::GetNumProcs()<3)
		{
			std::set<unsigned> local_boxes_to_box_1 = box_collection.GetLocalBoxes(1);

			std::set<unsigned> correct_answer_1;
			correct_answer_1.insert(1);
			correct_answer_1.insert(2);
			correct_answer_1.insert(3);
			correct_answer_1.insert(4);
			correct_answer_1.insert(5);
			correct_answer_1.insert(10);
			correct_answer_1.insert(11);
			correct_answer_1.insert(12);
			correct_answer_1.insert(13);
			correct_answer_1.insert(14);
			TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);
		}
		if(box_collection.GetBoxOwnership(1) && PetscTools::GetNumProcs()==3)
		{
			std::set<unsigned> local_boxes_to_box_1 = box_collection.GetLocalBoxes(1);

			std::set<unsigned> correct_answer_1;
			correct_answer_1.insert(0);
			correct_answer_1.insert(1);
			correct_answer_1.insert(2);
			correct_answer_1.insert(3);
			correct_answer_1.insert(4);
			correct_answer_1.insert(5);
			correct_answer_1.insert(9);
			correct_answer_1.insert(10);
			correct_answer_1.insert(11);
			correct_answer_1.insert(12);
			correct_answer_1.insert(13);
			correct_answer_1.insert(14);
			TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);
		}

		if(box_collection.GetBoxOwnership(2) && PetscTools::GetNumProcs()==1)
		{
			std::set<unsigned> local_boxes_to_box_2 = box_collection.GetLocalBoxes(2);

			std::set<unsigned> correct_answer_2;
			correct_answer_2.insert(2);
			correct_answer_2.insert(4);
			correct_answer_2.insert(5);
			correct_answer_2.insert(11);
			correct_answer_2.insert(14);
			correct_answer_2.insert(13);

			TS_ASSERT_EQUALS(local_boxes_to_box_2, correct_answer_2);
		}
		if(box_collection.GetBoxOwnership(2) && PetscTools::GetNumProcs()>1)
		{
			std::set<unsigned> local_boxes_to_box_2 = box_collection.GetLocalBoxes(2);

			std::set<unsigned> correct_answer_2;
			correct_answer_2.insert(1);
			correct_answer_2.insert(2);
			correct_answer_2.insert(4);
			correct_answer_2.insert(5);
			correct_answer_2.insert(10);
			correct_answer_2.insert(11);
			correct_answer_2.insert(14);
			correct_answer_2.insert(13);

			TS_ASSERT_EQUALS(local_boxes_to_box_2, correct_answer_2);
		}

		if(box_collection.GetBoxOwnership(17) && PetscTools::GetNumProcs()==1)
		{
			std::set<unsigned> local_boxes_to_box_17 = box_collection.GetLocalBoxes(17);

			std::set<unsigned> correct_answer_17;
			correct_answer_17.insert(8);
			correct_answer_17.insert(17);

			TS_ASSERT_EQUALS(local_boxes_to_box_17, correct_answer_17);
		}
		if(box_collection.GetBoxOwnership(17) && PetscTools::GetNumProcs()>1)
		{
			std::set<unsigned> local_boxes_to_box_17 = box_collection.GetLocalBoxes(17);

			std::set<unsigned> correct_answer_17;
			correct_answer_17.insert(17);
			correct_answer_17.insert(16);
			correct_answer_17.insert(13);
			correct_answer_17.insert(8);
			correct_answer_17.insert(7);
			correct_answer_17.insert(4);

			TS_ASSERT_EQUALS(local_boxes_to_box_17, correct_answer_17);
		}



	}

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    //
    //	Parallelisation Tests
    //	The following tests are written to test the parallel implementaiton of the
	//	box collection.
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    void TestSetupHaloBoxes1d2d3d() throw (Exception)
	{
    	// np == 1 - All halos should be empty
    	if(PetscTools::IsSequential())
    	{
			// 1D ///////////////////////////////////////////////////////////////////

			double cut_off_length = 5.0;
			c_vector<double, 2> domain_size1d;
			domain_size1d(0) = 0.0;
			domain_size1d(1) = 20.0;

			BoxCollection<1> box_collection1d(cut_off_length, domain_size1d);

			box_collection1d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right1d;
			std::vector<unsigned> halos_should_be_left1d;

			TS_ASSERT_EQUALS(halos_should_be_right1d, box_collection1d.mHalosRight);
			TS_ASSERT_EQUALS(halos_should_be_left1d, box_collection1d.mHalosLeft);
			TS_ASSERT_EQUALS(box_collection1d.mHaloBoxes.size(),0u);

			// 2D ///////////////////////////////////////////////////////////////////

			// Set up a box collection
			c_vector<double, 4> domain_size2d;
			domain_size2d(0) = 0.0;
			domain_size2d(1) = 20.0;
			domain_size2d(2) = 0.0;
			domain_size2d(3) = 20.0;

			BoxCollection<2> box_collection2d(cut_off_length, domain_size2d);

			box_collection2d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right2d;
			std::vector<unsigned> halos_should_be_left2d;

			TS_ASSERT_EQUALS(halos_should_be_right2d, box_collection2d.mHalosRight);
			TS_ASSERT_EQUALS(halos_should_be_left2d, box_collection2d.mHalosLeft);
			TS_ASSERT_EQUALS(box_collection2d.mHaloBoxes.size(),0u);

			// 3D ///////////////////////////////////////////////////////////////////

			// Set up a box collection
			c_vector<double, 6> domain_size3d;
			domain_size3d(0) = 0.0;
			domain_size3d(1) = 20.0;
			domain_size3d(2) = 0.0;
			domain_size3d(3) = 20.0;
			domain_size3d(4) = 0.0;
			domain_size3d(5) = 20.0;

			BoxCollection<3> box_collection3d(cut_off_length, domain_size3d);

			box_collection3d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right3d;
			std::vector<unsigned> halos_should_be_left3d;

			TS_ASSERT_EQUALS(halos_should_be_right3d, box_collection3d.mHalosRight);
			TS_ASSERT_EQUALS(halos_should_be_left3d, box_collection3d.mHalosLeft);
			TS_ASSERT_EQUALS(box_collection3d.mHaloBoxes.size(),0u);
    	}

    	if(PetscTools::GetNumProcs()==2)
		{
			// 1D ///////////////////////////////////////////////////////////////////

			double cut_off_length = 5.0;
			c_vector<double, 2> domain_size1d;
			domain_size1d(0) = 0.0;
			domain_size1d(1) = 20.0;

			BoxCollection<1> box_collection1d(cut_off_length, domain_size1d);

			box_collection1d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right1d;
			std::vector<unsigned> halos_should_be_left1d;



			if(PetscTools::GetMyRank()==0)
			{
				halos_should_be_right1d.push_back(1);
				TS_ASSERT_EQUALS(halos_should_be_right1d, box_collection1d.mHalosRight);
			}
			if(PetscTools::GetMyRank()==1)
			{
				halos_should_be_left1d.push_back(2);
				TS_ASSERT_EQUALS(halos_should_be_left1d, box_collection1d.mHalosLeft);
			}

			// 2D ///////////////////////////////////////////////////////////////////

			// Set up a box collection
			c_vector<double, 4> domain_size2d;
			domain_size2d(0) = 0.0;
			domain_size2d(1) = 20.0;
			domain_size2d(2) = 0.0;
			domain_size2d(3) = 20.0;

			BoxCollection<2> box_collection2d(cut_off_length, domain_size2d);

			box_collection2d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right2d;
			std::vector<unsigned> halos_should_be_left2d;

			if(PetscTools::GetMyRank()==0)
			{
				for(unsigned i=0;i<4;i++)
				{
					halos_should_be_right2d.push_back(1+4*i+2*PetscTools::GetMyRank());
				}
				TS_ASSERT_EQUALS(halos_should_be_right2d, box_collection2d.mHalosRight);
			}
			if(PetscTools::GetMyRank()==1)
			{
				for(unsigned i=0;i<4;i++)
				{
					halos_should_be_left2d.push_back(4*i+2*PetscTools::GetMyRank());
				}
				TS_ASSERT_EQUALS(halos_should_be_left2d, box_collection2d.mHalosLeft);
			}

			// 3D ///////////////////////////////////////////////////////////////////

			// Set up a box collection
			c_vector<double, 6> domain_size3d;
			domain_size3d(0) = 0.0;
			domain_size3d(1) = 20.0;
			domain_size3d(2) = 0.0;
			domain_size3d(3) = 20.0;
			domain_size3d(4) = 0.0;
			domain_size3d(5) = 20.0;

			BoxCollection<3> box_collection3d(cut_off_length, domain_size3d);

			box_collection3d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right3d;
			std::vector<unsigned> halos_should_be_left3d;

			if(PetscTools::GetMyRank() == 0)
			{
				for(unsigned j=0;j<4;j++)
				{
					for(unsigned i=0;i<4;i++)
					{
						halos_should_be_right3d.push_back((1+4*i)+16*j+2*PetscTools::GetMyRank());
					}
				}
				TS_ASSERT_EQUALS(halos_should_be_right3d, box_collection3d.mHalosRight);
			}
			if(PetscTools::GetMyRank() == 1)
			{
				for(unsigned j=0;j<4;j++)
				{
					for(unsigned i=0;i<4;i++)
					{
						halos_should_be_left3d.push_back((4*i)+16*j+2*PetscTools::GetMyRank());
					}
				}
				TS_ASSERT_EQUALS(halos_should_be_left3d, box_collection3d.mHalosLeft);
			}


		}

    	if(PetscTools::GetNumProcs()==3)
		{
			// 1D ///////////////////////////////////////////////////////////////////

			double cut_off_length = 5.0;
			c_vector<double, 2> domain_size1d;
			domain_size1d(0) = 0.0;
			domain_size1d(1) = 15.0;

			BoxCollection<1> box_collection1d(cut_off_length, domain_size1d);

			box_collection1d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right1d;
			std::vector<unsigned> halos_should_be_left1d;

			if(!PetscTools::AmMaster())
			{
				halos_should_be_left1d.push_back(PetscTools::GetMyRank());
			}
			if(!PetscTools::AmTopMost())
			{
				halos_should_be_right1d.push_back(PetscTools::GetMyRank());
			}

			TS_ASSERT_EQUALS(halos_should_be_left1d, box_collection1d.mHalosLeft);
			TS_ASSERT_EQUALS(halos_should_be_right1d, box_collection1d.mHalosRight);


			// 2D ///////////////////////////////////////////////////////////////////

			// Set up a box collection
			c_vector<double, 4> domain_size2d;
			domain_size2d(0) = 0.0;
			domain_size2d(1) = 15.0;
			domain_size2d(2) = 0.0;
			domain_size2d(3) = 15.0;

			BoxCollection<2> box_collection2d(cut_off_length, domain_size2d);

			box_collection2d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right2d;
			std::vector<unsigned> halos_should_be_left2d;

			if(PetscTools::GetMyRank()<2)
			{
				for(unsigned i=0;i<3;i++)
				{
					halos_should_be_right2d.push_back(PetscTools::GetMyRank() + 3*i);
				}
			}

			if(PetscTools::GetMyRank()>0)
			{
				for(unsigned i=0;i<3;i++)
				{
					halos_should_be_left2d.push_back(PetscTools::GetMyRank() + 3*i);
				}
			}
			TS_ASSERT_EQUALS(halos_should_be_right2d, box_collection2d.mHalosRight);
			TS_ASSERT_EQUALS(halos_should_be_left2d, box_collection2d.mHalosLeft);

			// 3D ///////////////////////////////////////////////////////////////////

			// Set up a box collection
			c_vector<double, 6> domain_size3d;
			domain_size3d(0) = 0.0;
			domain_size3d(1) = 20.0;
			domain_size3d(2) = 0.0;
			domain_size3d(3) = 20.0;
			domain_size3d(4) = 0.0;
			domain_size3d(5) = 20.0;

			BoxCollection<3> box_collection3d(cut_off_length, domain_size3d);

			box_collection3d.SetupHaloBoxes();

			std::vector<unsigned> halos_should_be_right3d;
			std::vector<unsigned> halos_should_be_left3d;
			if(PetscTools::GetMyRank()<2)
			{
				for(unsigned j=0;j<4;j++)
				{
					for(unsigned i=0;i<4;i++)
					{
						halos_should_be_right3d.push_back(PetscTools::GetMyRank()+1+4*i+16*j);
					}
				}
			}
			if(PetscTools::GetMyRank()>0)
			{
				for(unsigned j=0;j<4;j++)
				{
					for(unsigned i=0;i<4;i++)
					{
						halos_should_be_left3d.push_back(PetscTools::GetMyRank()+1 + (4*i)+16*j);
					}
				}
			}

			TS_ASSERT_EQUALS(halos_should_be_right3d, box_collection3d.mHalosRight);
			TS_ASSERT_EQUALS(halos_should_be_left3d, box_collection3d.mHalosLeft);
		}
	}

    void TestUpdateHaloBoxes1d2d3d() throw (Exception)
	{
    	if(PetscTools::GetNumProcs()==2)
    	{
    		//////////////////
    		// 1 Dimension	//
    		//////////////////
    		{
				std::vector< ChastePoint<1>* > points(2);
				points[0] = new ChastePoint<1>(0.5);
				points[1] = new ChastePoint<1>(1.5);

				std::vector<Node<1>* > nodes;
				for (unsigned i=0; i<points.size(); i++)
				{
					nodes.push_back(new Node<1>(i, *(points[i]), false));
				}

				double cut_off_length = 1.0;

				c_vector<double, 2> domain_size;
				domain_size(0) = 0.0;
				domain_size(1) = 2.0;

				BoxCollection<1> box_collection(cut_off_length, domain_size);
				box_collection.SetupHaloBoxes();

				if(PetscTools::GetMyRank()==0)
				{
					box_collection.mBoxes[0].AddNode(nodes[0]);
				}
				if(PetscTools::GetMyRank()==1)
				{
					box_collection.mBoxes[0].AddNode(nodes[1]);
				}

				box_collection.UpdateHaloBoxes();
				// Coverage of deleting nodes when updating.
				box_collection.UpdateHaloBoxes();

				// Check correct nodes lie in halos.
				for(unsigned i=0;i<box_collection.mHaloBoxes.size();i++)
				{
					std::set< Node<1>* > halonodes = box_collection.mHaloBoxes[i].rGetNodesContained();

					// Iterate over the nodes
					for( std::set< Node<1>* >::iterator it=halonodes.begin();
							it!=halonodes.end();
							it++)
					{
						c_vector<double, 1> location = (*it)->rGetLocation();
						c_vector<double, 1> correct_location;
						correct_location[0]=1.5-PetscTools::GetMyRank();
						TS_ASSERT_EQUALS(location[0],correct_location[0]);
					}
				}
    		}
    		//////////////////
			// 2 Dimensions	//
			//////////////////
    		{
				std::vector< ChastePoint<2>* > points(4);
				points[0] = new ChastePoint<2>(0.5,0.5);
				points[1] = new ChastePoint<2>(1.5,0.5);
				points[2] = new ChastePoint<2>(0.5,1.5);
				points[3] = new ChastePoint<2>(1.5,1.5);

				std::vector<Node<2>* > nodes;
				for (unsigned i=0; i<points.size(); i++)
				{
					nodes.push_back(new Node<2>(i, *(points[i]), false));
				}

				double cut_off_length = 1.0;

				c_vector<double, 4> domain_size;
				domain_size(0) = 0.0;
				domain_size(1) = 2.0;
				domain_size(2) = 0.0;
				domain_size(3) = 2.0;

				BoxCollection<2> box_collection(cut_off_length, domain_size);
				box_collection.SetupHaloBoxes();

				if(PetscTools::GetMyRank()==0)
				{
					box_collection.mBoxes[0].AddNode(nodes[0]);
					box_collection.mBoxes[1].AddNode(nodes[2]);
				}
				if(PetscTools::GetMyRank()==1)
				{
					box_collection.mBoxes[0].AddNode(nodes[1]);
					box_collection.mBoxes[1].AddNode(nodes[3]);
				}

				box_collection.UpdateHaloBoxes();

				// Check correct nodes lie in halos.
				for(unsigned i=0;i<box_collection.mHaloBoxes.size();i++)
				{
					std::set< Node<2>* > halonodes = box_collection.mHaloBoxes[i].rGetNodesContained();

					// Iterate over the nodes
					for( std::set< Node<2>* >::iterator it=halonodes.begin();
							it!=halonodes.end();
							it++)
					{
						c_vector<double, 2> location = (*it)->rGetLocation();
						c_vector<double, 2> correct_location;
						correct_location[0]=1.5-PetscTools::GetMyRank();
						correct_location[1]=0.5+(double)i;
						TS_ASSERT_EQUALS(location[0],correct_location[0]);
						TS_ASSERT_EQUALS(location[1],correct_location[1]);
					}
				}
    		}
    		//////////////////
			// 3 Dimensions	//
			//////////////////
			{
				std::vector< ChastePoint<3>* > points(8);
				points[0] = new ChastePoint<3>(0.5,0.5,0.5);
				points[1] = new ChastePoint<3>(1.5,0.5,0.5);
				points[2] = new ChastePoint<3>(0.5,1.5,0.5);
				points[3] = new ChastePoint<3>(1.5,1.5,0.5);
				points[4] = new ChastePoint<3>(0.5,0.5,1.5);
				points[5] = new ChastePoint<3>(1.5,0.5,1.5);
				points[6] = new ChastePoint<3>(0.5,1.5,1.5);
				points[7] = new ChastePoint<3>(1.5,1.5,1.5);

				std::vector<Node<3>* > nodes;
				for (unsigned i=0; i<points.size(); i++)
				{
					nodes.push_back(new Node<3>(i, *(points[i]), false));
				}

				double cut_off_length = 1.0;

				c_vector<double, 6> domain_size;
				domain_size(0) = 0.0;
				domain_size(1) = 2.0;
				domain_size(2) = 0.0;
				domain_size(3) = 2.0;
				domain_size(4) = 0.0;
				domain_size(5) = 2.0;

				BoxCollection<3> box_collection(cut_off_length, domain_size);
				box_collection.SetupHaloBoxes();

				if(PetscTools::GetMyRank()==0)
				{
					box_collection.mBoxes[0].AddNode(nodes[0]);
					box_collection.mBoxes[1].AddNode(nodes[2]);
					box_collection.mBoxes[2].AddNode(nodes[4]);
					box_collection.mBoxes[3].AddNode(nodes[6]);
				}
				if(PetscTools::GetMyRank()==1)
				{
					box_collection.mBoxes[0].AddNode(nodes[1]);
					box_collection.mBoxes[1].AddNode(nodes[3]);
					box_collection.mBoxes[2].AddNode(nodes[5]);
					box_collection.mBoxes[3].AddNode(nodes[7]);
				}

				box_collection.UpdateHaloBoxes();

				// Check correct nodes lie in halos.
				for(unsigned i=0;i<box_collection.mHaloBoxes.size();i++)
				{
					std::set< Node<3>* > halonodes = box_collection.mHaloBoxes[i].rGetNodesContained();

					// Iterate over the nodes
					for( std::set< Node<3>* >::iterator it=halonodes.begin();
							it!=halonodes.end();
							it++)
					{
						c_vector<double, 3> location = (*it)->rGetLocation();
						c_vector<double, 3> correct_location;

						correct_location[0]=1.5-PetscTools::GetMyRank();
						correct_location[1]=0.5+(double)(i%2);
						correct_location[2]=0.5+(double)(unsigned)(i/2);
						TS_ASSERT_EQUALS(location[0],correct_location[0]);
						TS_ASSERT_EQUALS(location[1],correct_location[1]);
						TS_ASSERT_EQUALS(location[2],correct_location[2]);
					}
				}
			}
    	}
		if(PetscTools::GetNumProcs()==3)
		{
			//////////////////
			// 1 Dimension	//
			//////////////////
			{
				std::vector< ChastePoint<1>* > points(3);
				points[0] = new ChastePoint<1>(0.5);
				points[1] = new ChastePoint<1>(1.5);
				points[2] = new ChastePoint<1>(2.5);

				std::vector<Node<1>* > nodes;
				for (unsigned i=0; i<points.size(); i++)
				{
					nodes.push_back(new Node<1>(i, *(points[i]), false));
				}

				double cut_off_length = 1.0;

				c_vector<double, 2> domain_size;
				domain_size(0) = 0.0;
				domain_size(1) = 3.0;

				BoxCollection<1> box_collection(cut_off_length, domain_size);
				box_collection.SetupHaloBoxes();

				if(PetscTools::GetMyRank()==0)
				{
					box_collection.mBoxes[0].AddNode(nodes[0]);
				}
				if(PetscTools::GetMyRank()==1)
				{
					box_collection.mBoxes[0].AddNode(nodes[1]);
				}
				if(PetscTools::GetMyRank()==2)
				{
					box_collection.mBoxes[0].AddNode(nodes[2]);
				}

				box_collection.UpdateHaloBoxes();

				// Check correct nodes lie in halos.
				for(unsigned i=0;i<box_collection.mHaloBoxes.size();i++)
				{
					std::set< Node<1>* > halonodes = box_collection.mHaloBoxes[i].rGetNodesContained();

					// Iterate over the nodes
					for( std::set< Node<1>* >::iterator it=halonodes.begin();
							it!=halonodes.end();
							it++)
					{
						if(PetscTools::AmMaster())
						{
							c_vector<double, 1> location = (*it)->rGetLocation();
							c_vector<double, 1> correct_location;
							correct_location[0]=1.5;
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
						}
						else if(PetscTools::AmTopMost())
						{
							c_vector<double, 1> location = (*it)->rGetLocation();
							c_vector<double, 1> correct_location;
							correct_location[0]=1.5;
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
						}
						else if(PetscTools::AmTopMost())
						{
							c_vector<double, 1> location = (*it)->rGetLocation();
							c_vector<double, 1> correct_location;
							correct_location[0]=0.5+2*(unsigned)i;
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
						}
					}
				}
			}
			//////////////////
			// 2 Dimensions	//
			//////////////////
			{
				std::vector< ChastePoint<2>* > points(9);
				points[0] = new ChastePoint<2>(0.5,0.5);
				points[1] = new ChastePoint<2>(1.5,0.5);
				points[2] = new ChastePoint<2>(2.5,0.5);
				points[3] = new ChastePoint<2>(0.5,1.5);
				points[4] = new ChastePoint<2>(1.5,1.5);
				points[5] = new ChastePoint<2>(2.5,1.5);
				points[6] = new ChastePoint<2>(0.5,2.5);
				points[7] = new ChastePoint<2>(1.5,2.5);
				points[8] = new ChastePoint<2>(2.5,2.5);

				std::vector<Node<2>* > nodes;
				for (unsigned i=0; i<points.size(); i++)
				{
					nodes.push_back(new Node<2>(i, *(points[i]), false));
				}

				double cut_off_length = 1.0;

				c_vector<double, 4> domain_size;
				domain_size(0) = 0.0;
				domain_size(1) = 3.0;
				domain_size(2) = 0.0;
				domain_size(3) = 3.0;

				BoxCollection<2> box_collection(cut_off_length, domain_size);
				box_collection.SetupHaloBoxes();

				if(PetscTools::GetMyRank()==0)
				{
					box_collection.mBoxes[0].AddNode(nodes[0]);
					box_collection.mBoxes[1].AddNode(nodes[3]);
					box_collection.mBoxes[2].AddNode(nodes[6]);
				}
				if(PetscTools::GetMyRank()==1)
				{
					box_collection.mBoxes[0].AddNode(nodes[1]);
					box_collection.mBoxes[1].AddNode(nodes[4]);
					box_collection.mBoxes[2].AddNode(nodes[7]);
				}
				if(PetscTools::GetMyRank()==2)
				{
					box_collection.mBoxes[0].AddNode(nodes[2]);
					box_collection.mBoxes[1].AddNode(nodes[5]);
					box_collection.mBoxes[2].AddNode(nodes[8]);
				}

				box_collection.UpdateHaloBoxes();

				// Check correct nodes lie in halos.
				for(unsigned i=0;i<box_collection.mHaloBoxes.size();i++)
				{
					std::set< Node<2>* > halonodes = box_collection.mHaloBoxes[i].rGetNodesContained();

					// Iterate over the nodes
					for( std::set< Node<2>* >::iterator it=halonodes.begin();
							it!=halonodes.end();
							it++)
					{
						if(PetscTools::AmMaster())
						{
							c_vector<double, 2> location = (*it)->rGetLocation();
							c_vector<double, 2> correct_location;
							correct_location[0]=1.5;
							correct_location[1]=0.5+(double)i;
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
							TS_ASSERT_EQUALS(location[1],correct_location[1]);
						}
						else if(PetscTools::AmTopMost())
						{
							c_vector<double, 2> location = (*it)->rGetLocation();
							c_vector<double, 2> correct_location;
							correct_location[0]=1.5;
							correct_location[1]=0.5+(double)i;
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
							TS_ASSERT_EQUALS(location[1],correct_location[1]);
						}
						else
						{
							c_vector<double, 2> location = (*it)->rGetLocation();
							c_vector<double, 2> correct_location;
							if(i<3)
							{
								correct_location[0]=2.5;
								correct_location[1]=0.5+(double)i;
							}
							else
							{
								correct_location[0]=0.5;
								correct_location[1]=0.5+(double)(i-3);
							}

							TS_ASSERT_EQUALS(location[0],correct_location[0]);
							TS_ASSERT_EQUALS(location[1],correct_location[1]);
						}
					}
				}
			}
			//////////////////
			// 3 Dimensions	//
			//////////////////
			{
				std::vector< ChastePoint<3>* > points(27);
				for(unsigned i=0; i<3; i++)
				{
					double x = 0.5+(double)i;

					for(unsigned j=0; j<3; j++)
					{
						double y = 0.5+ (double)j;

						for(unsigned k=0; k<3; k++)
						{
							double z = 0.5+ (double)k;
							points[k+3*j+9*i] = new ChastePoint<3>(x,y,z);
						}
					}
				}

				std::vector<Node<3>* > nodes;
				for (unsigned i=0; i<points.size(); i++)
				{
					nodes.push_back(new Node<3>(i, *(points[i]), false));
				}

				double cut_off_length = 1.0;

				c_vector<double, 6> domain_size;
				domain_size(0) = 0.0;
				domain_size(1) = 3.0;
				domain_size(2) = 0.0;
				domain_size(3) = 3.0;
				domain_size(4) = 0.0;
				domain_size(5) = 3.0;

				BoxCollection<3> box_collection(cut_off_length, domain_size);
				box_collection.SetupHaloBoxes();

				if(PetscTools::GetMyRank()==0)
				{
					for(unsigned i =0; i < 9; i++)
					{
						box_collection.mBoxes[i].AddNode(nodes[i]);
					}
				}
				if(PetscTools::GetMyRank()==1)
				{
					for(unsigned i =0; i < 9; i++)
					{
						box_collection.mBoxes[i].AddNode(nodes[i+9]);
					}
				}
				if(PetscTools::GetMyRank()==2)
				{
					for(unsigned i =0; i < 9; i++)
					{
						box_collection.mBoxes[i].AddNode(nodes[i+18]);
					}
				}

				box_collection.UpdateHaloBoxes();

				// Check correct nodes lie in halos.
				for(unsigned i=0;i<box_collection.mHaloBoxes.size();i++)
				{
					std::set< Node<3>* > halonodes = box_collection.mHaloBoxes[i].rGetNodesContained();

					// Iterate over the nodes
					for( std::set< Node<3>* >::iterator it=halonodes.begin();
							it!=halonodes.end();
							it++)
					{
						if(PetscTools::AmMaster())
						{
							c_vector<double, 3> location = (*it)->rGetLocation();
							c_vector<double, 3> correct_location;
							correct_location[0]=1.5;
							correct_location[1]=0.5+(double)(i%3);
							correct_location[2]=0.5+(double)(unsigned)(i/3);
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
							TS_ASSERT_EQUALS(location[1],correct_location[1]);
							TS_ASSERT_EQUALS(location[2],correct_location[2]);
						}
						else if(PetscTools::AmTopMost())
						{
							c_vector<double, 3> location = (*it)->rGetLocation();
							c_vector<double, 3> correct_location;
							correct_location[0]=1.5;
							correct_location[1]=0.5+(double)(i%3);
							correct_location[2]=0.5+(double)(unsigned)(i/3);
							TS_ASSERT_EQUALS(location[0],correct_location[0]);
							TS_ASSERT_EQUALS(location[1],correct_location[1]);
							TS_ASSERT_EQUALS(location[2],correct_location[2]);
						}
						else
						{
							c_vector<double, 3> location = (*it)->rGetLocation();
							c_vector<double, 3> correct_location;
							if(i<9)
							{
								correct_location[0]=2.5;
								correct_location[1]=0.5+(double)(i%3);
								correct_location[2]=0.5+(double)(unsigned)(i/3);
							}
							else
							{
								correct_location[0]=0.5;
								correct_location[1]=0.5+(double)((i-3)%3);
								correct_location[2]=0.5+(double)(unsigned)((i-9)/3);
							}

							TS_ASSERT_EQUALS(location[0],correct_location[0]);
							TS_ASSERT_EQUALS(location[1],correct_location[1]);
							TS_ASSERT_EQUALS(location[2],correct_location[2]);
						}
					}
				}
			}
		}
	}


    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    //
    //  Cancer/cell_based tests
    //  The following are tests written from this BoxCollection used to be
    //  cell_based/src/tissue/NodeBoxCollection and test the cell-based
    //  functionality
    //
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    void TestPairsReturned1d() throw (Exception)
    {
		std::vector< ChastePoint<1>* > points(7);
		points[0] = new ChastePoint<1>(0.5);
		points[1] = new ChastePoint<1>(1.5);
		points[2] = new ChastePoint<1>(2.5);
		points[3] = new ChastePoint<1>(3.5);
		points[4] = new ChastePoint<1>(4.5);
		points[5] = new ChastePoint<1>(5.5);
		points[6] = new ChastePoint<1>(5.6);


		std::vector<Node<1>* > nodes;
		for (unsigned i=0; i<points.size(); i++)
		{
			nodes.push_back(new Node<1>(i, *(points[i]), false));
		}

		double cut_off_length = 1.0;

		c_vector<double, 2> domain_size;
		domain_size(0) = 0.0;
		domain_size(1) = 6.0;

		BoxCollection<1> box_collection(cut_off_length, domain_size);

		box_collection.SetupLocalBoxesHalfOnly();
		TS_ASSERT_THROWS_THIS(box_collection.SetupLocalBoxesHalfOnly(), "Local Boxes Are Already Set");

		box_collection.SetupHaloBoxes();

		for (unsigned i=0; i<nodes.size(); i++)
		{
			unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
			if(box_collection.GetBoxOwnership(box_index))
			{
				box_collection.rGetBox(box_index).AddNode(nodes[i]);
			}
		}


		std::set< std::pair<Node<1>*, Node<1>* > > pairs_returned;
		box_collection.CalculateNodePairs(nodes,pairs_returned);

		std::set< std::pair<Node<1>*, Node<1>* > > pairs_should_be;

		// Set out what answers should be
		if(PetscTools::IsSequential())
		{
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[3],nodes[4]));
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[5]));
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[5],nodes[6]));
			pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[6]));

			TS_ASSERT_EQUALS(box_collection.SolveBoxMapping(5),5u);

		}
		else if(PetscTools::GetNumProcs()==2)
		{

			if(PetscTools::GetMyRank()==0)
			{
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));
			}
			else if(PetscTools::GetMyRank()==1)
			{
				// Coverage
				TS_ASSERT_THROWS_THIS(box_collection.rGetBox(0), "Box does not exist on this process! Try calling GetBoxOwnership() before rGetBox.");
				TS_ASSERT_EQUALS(box_collection.SolveBoxMapping(5),2u);

				TS_ASSERT_THROWS_THIS(box_collection.SolveBoxMapping(0), "Requested box with global index 0, which does not belong to processor 1");

				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[3],nodes[2]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[3],nodes[4]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[5]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[6]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[5],nodes[6]));
			}
		}
		else if(PetscTools::GetNumProcs()==3)
		{
			if(PetscTools::GetMyRank()==0)
			{
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
			}
			else if(PetscTools::GetMyRank()==1)
			{
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[1]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[3],nodes[4]));
			}
			else if(PetscTools::GetMyRank()==2)
			{
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[3]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[5]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[4],nodes[6]));
				pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[5],nodes[6]));
				TS_ASSERT_EQUALS(box_collection.SolveBoxMapping(5),1u);
			}


		}

		TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

		for (unsigned i=0; i<points.size(); i++)
		{
			delete nodes[i];
			delete points[i];
		}
    }

    void TestBoxGeneration2d() throw (Exception)
    {
    	// Create a simple mesh
		TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
		TetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

		double cut_off_length = 0.2;

		c_vector<double, 2*2> domain_size;
		domain_size(0) = -0.1;
		domain_size(1) = 1.15;
		domain_size(2) = -0.1;
		domain_size(3) = 1.15;

		BoxCollection<2> box_collection(cut_off_length, domain_size);

		box_collection.SetupLocalBoxesHalfOnly();

		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
			if(box_collection.GetBoxOwnership(box_index))
			{
				box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
			}
		}

		TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 49u);

		for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
		{
			if(box_collection.GetBoxOwnership(i))
			{
				std::set< Node<2>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
				c_vector<double, 2*2> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

				for (std::set< Node<2>* >::iterator it_nodes_in_box = nodes_in_box.begin();
					 it_nodes_in_box != nodes_in_box.end();
					 it_nodes_in_box++)
				{
					Node<2>* current_node = *it_nodes_in_box;
					double x_position = current_node->rGetLocation()[0];
					double y_position = current_node->rGetLocation()[1];

					double epsilon = 1e-12;

					TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
					TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
					TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon, y_position);
					TS_ASSERT_LESS_THAN(y_position, box_min_max_values(3)+epsilon);
				}
			}
		}

		// Have checked that all the local boxes are calculated correctly on a 5 by 6 grid - here we
		// hardcode a few checks on the 7 by 7 grid.
		if(PetscTools::IsSequential())
		{
			std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
			std::set<unsigned> correct_answer_0;
			correct_answer_0.insert(0);
			correct_answer_0.insert(1);
			correct_answer_0.insert(7);
			correct_answer_0.insert(8);
			TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

			std::set<unsigned> local_boxes_to_box_4 = box_collection.GetLocalBoxes(4);
			std::set<unsigned> correct_answer_4;
			correct_answer_4.insert(4);
			correct_answer_4.insert(5);
			correct_answer_4.insert(10);
			correct_answer_4.insert(11);
			correct_answer_4.insert(12);
			TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

			std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
			std::set<unsigned> correct_answer_10;
			correct_answer_10.insert(10);
			correct_answer_10.insert(11);
			correct_answer_10.insert(16);
			correct_answer_10.insert(17);
			correct_answer_10.insert(18);
			TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

			std::set<unsigned> local_boxes_to_box_48 = box_collection.GetLocalBoxes(48);
			std::set<unsigned> correct_answer_48;
			correct_answer_48.insert(48);
			TS_ASSERT_EQUALS(local_boxes_to_box_48, correct_answer_48);
		}

		c_vector<unsigned, 2> indices=box_collection.CalculateCoordinateIndices(0);
		TS_ASSERT_EQUALS(indices[0],0u);
		TS_ASSERT_EQUALS(indices[1],0u);
    }

    void TestPairsReturned2d() throw (Exception)
    {
		std::vector< ChastePoint<2>* > points(6);

		points[0]= new ChastePoint<2>(0.5, 0.5);
		points[1]= new ChastePoint<2>(1.5, 0.5);
		points[2]= new ChastePoint<2>(0.5, 1.5);
		points[3]= new ChastePoint<2>(2.5, 1.5);
		points[4]= new ChastePoint<2>(1.5, 2.5);
		points[5]= new ChastePoint<2>(2.5, 2.5);

		std::vector<Node<2>* > nodes;
		for (unsigned i=0; i<points.size(); i++)
		{
			nodes.push_back(new Node<2>(i, *(points[i]), false));
		}

		double cut_off_length = 1.0;

		c_vector<double, 2*2> domain_size;
		domain_size(0) = 0.0;
		domain_size(1) = 3.0;
		domain_size(2) = 0.0;
		domain_size(3) = 3.0;

		BoxCollection<2> box_collection(cut_off_length, domain_size);

		box_collection.SetupLocalBoxesHalfOnly();
		box_collection.SetupHaloBoxes();

		std::set< std::pair<Node<2>*, Node<2>* > > pairs_returned;
		box_collection.CalculateNodePairs(nodes,pairs_returned);


		std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
		if(PetscTools::GetNumProcs()==1)
		{
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[2]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[3]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[2]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[4]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
			pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
		}
		if(PetscTools::GetNumProcs()==2)
		{
			if(PetscTools::GetMyRank()==0)
			{
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[2]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[2]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[3]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[4]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[3]));
			}
			if(PetscTools::GetMyRank()==1)
			{
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[1]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[4]));
			}
		}
		if(PetscTools::GetNumProcs()==3)
		{
			if(PetscTools::GetMyRank()==0)
			{
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[2]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[1]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[4]));
			}

			if(PetscTools::GetMyRank()==1)
			{
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[0]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[2]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[1],nodes[3]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[3]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[2]));
			}
			if(PetscTools::GetMyRank()==2)
			{
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[1]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
				pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[4]));
			}
		}

		TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

		for (unsigned i=0; i<points.size(); i++)
		{
			delete nodes[i];
			delete points[i];
		}
    }


    void TestBoxGeneration3d() throw (Exception)
    {
    	//Currently on tested on 1 process.
    	if(PetscTools::IsSequential())
    	{
			// Create a mesh
			TetrahedralMesh<3,3> mesh;
			mesh.ConstructCuboid(4,5,6);

			double cut_off_length = 2.0;

			c_vector<double, 2*3> domain_size;
			domain_size(0) = -0.1;
			domain_size(1) = 4.15;
			domain_size(2) = -0.1;
			domain_size(3) = 5.15;
			domain_size(4) = -0.1;
			domain_size(5) = 6.15;

			BoxCollection<3> box_collection(cut_off_length, domain_size);

			box_collection.SetupLocalBoxesHalfOnly();

			for (unsigned i=0; i<mesh.GetNumNodes(); i++)
			{
				unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
				box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
			}

			TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 36u);

			for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
			{
				std::set< Node<3>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
				c_vector<double, 2*3> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

				for (std::set< Node<3>* >::iterator it_nodes_in_box = nodes_in_box.begin();
					 it_nodes_in_box != nodes_in_box.end();
					 it_nodes_in_box++)
				{
					Node<3>* current_node = *it_nodes_in_box;
					double x_position = current_node->rGetLocation()[0];
					double y_position = current_node->rGetLocation()[1];
					double z_position = current_node->rGetLocation()[2];

					double epsilon = 1e-12;

					TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
					TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
					TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon, y_position);
					TS_ASSERT_LESS_THAN(y_position, box_min_max_values(3)+epsilon);
					TS_ASSERT_LESS_THAN(box_min_max_values(4)-epsilon, z_position);
					TS_ASSERT_LESS_THAN(z_position, box_min_max_values(5)+epsilon);
				}
			}

			std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
			std::set<unsigned> correct_answer_0;
			correct_answer_0.insert(0);
			correct_answer_0.insert(1);
			correct_answer_0.insert(3);
			correct_answer_0.insert(4);
			correct_answer_0.insert(9);
			correct_answer_0.insert(10);
			correct_answer_0.insert(12);
			correct_answer_0.insert(13);
			TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

			std::set<unsigned> local_boxes_to_box_13 = box_collection.GetLocalBoxes(13);
			std::set<unsigned> correct_answer_13;
			correct_answer_13.insert(4);
			correct_answer_13.insert(5);
			correct_answer_13.insert(6);
			correct_answer_13.insert(7);
			correct_answer_13.insert(8);
			correct_answer_13.insert(13);
			correct_answer_13.insert(14);
			correct_answer_13.insert(15);
			correct_answer_13.insert(16);
			correct_answer_13.insert(17);
			correct_answer_13.insert(22);
			correct_answer_13.insert(23);
			correct_answer_13.insert(24);
			correct_answer_13.insert(25);
			correct_answer_13.insert(26);
			TS_ASSERT_EQUALS(local_boxes_to_box_13, correct_answer_13);

			std::set<unsigned> local_boxes_to_box_34 = box_collection.GetLocalBoxes(34);
			std::set<unsigned> correct_answer_34;
			correct_answer_34.insert(25);
			correct_answer_34.insert(26);
			correct_answer_34.insert(34);
			correct_answer_34.insert(35);
			TS_ASSERT_EQUALS(local_boxes_to_box_34, correct_answer_34);

			std::set<unsigned> local_boxes_to_box_35 = box_collection.GetLocalBoxes(35);
			std::set<unsigned> correct_answer_35;
			correct_answer_35.insert(26);
			correct_answer_35.insert(35);
			TS_ASSERT_EQUALS(local_boxes_to_box_35, correct_answer_35);


			c_vector<unsigned, 3> indices=box_collection.CalculateCoordinateIndices(0);
			TS_ASSERT_EQUALS(indices[0],0u);
			TS_ASSERT_EQUALS(indices[1],0u);
			TS_ASSERT_EQUALS(indices[2],0u);
		}

    }

    void TestPairsReturned3d() throw(Exception)
	{
    	std::vector< ChastePoint<3>* > points(9);

		points[0]= new ChastePoint<3>(0.5, 0.5, 0.5);
		points[1]= new ChastePoint<3>(1.5, 0.5, 0.5);
		points[2]= new ChastePoint<3>(0.5, 1.5, 0.5);
		points[3]= new ChastePoint<3>(2.5, 1.5, 0.5);
		points[4]= new ChastePoint<3>(1.5, 2.5, 0.5);
		points[5]= new ChastePoint<3>(2.5, 2.5, 0.5);
		points[6]= new ChastePoint<3>(2.5, 0.5, 1.5);
		points[7]= new ChastePoint<3>(1.5, 1.5, 1.5);
		points[8]= new ChastePoint<3>(0.5, 2.5, 1.5);

		std::vector<Node<3>* > nodes;
		for (unsigned i=0; i<points.size(); i++)
		{
			nodes.push_back(new Node<3>(i, *(points[i]), false));
		}

		double cut_off_length = 1.0;

		c_vector<double, 6> domain_size;
		domain_size(0) = 0.0;
		domain_size(1) = 3.0;
		domain_size(2) = 0.0;
		domain_size(3) = 3.0;
		domain_size(4) = 0.0;
		domain_size(5) = 2.0;

		BoxCollection<3> box_collection(cut_off_length, domain_size);

		box_collection.SetupLocalBoxesHalfOnly();
		box_collection.SetupHaloBoxes();

		std::set< std::pair<Node<3>*, Node<3>* > > pairs_returned;
		box_collection.CalculateNodePairs(nodes,pairs_returned);

		std::set< std::pair<Node<3>*, Node<3>* > > pairs_should_be;
		if(PetscTools::IsSequential())
		{
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[1]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[2]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[7]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[2]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[3]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[6]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[7]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[4]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[7]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[8]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[4]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[5]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[5]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[3]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[7]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[3]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[5]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[4]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[8]));
			pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[8],nodes[4]));
		}
		if(PetscTools::GetNumProcs()==2)
		{
			if(PetscTools::GetMyRank()==0)
			{
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[2]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[2]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[6]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[8]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[5]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[8]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[6]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[5]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[8],nodes[4]));
			}
			if(PetscTools::GetMyRank()==1)
			{
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[5]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[5],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[5],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[3]));

				// Coverage
				TS_ASSERT_THROWS_THIS(box_collection.GetLocalBoxes(0), "Called GetLocalBoxes on a boxIndex that does not belong to this process");
			}
		}
		if(PetscTools::GetNumProcs()==3)
		{
			if(PetscTools::GetMyRank()==0)
			{
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[2]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[0],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[2],nodes[8]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[8],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[8],nodes[7]));

			}
			if(PetscTools::GetMyRank()==1)
			{
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[0]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[2]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[6]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[1],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[2]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[5]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[4],nodes[8]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[0]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[2]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[5]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[6]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[7],nodes[8]));

			}
			if(PetscTools::GetMyRank()==2)
			{
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[5]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[3],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[5],nodes[4]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[5],nodes[7]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[1]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[3]));
				pairs_should_be.insert(std::pair<Node<3>*, Node<3>*>(nodes[6],nodes[7]));

			}
		}
		TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

		for (unsigned i=0; i<points.size(); i++)
		{
			delete nodes[i];
			delete points[i];
		}
	}
};

#endif /*TESTBOXCOLLECTION_HPP_*/
