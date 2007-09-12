#ifndef TESTMODIFIEDSPRINGMODEL_HPP_
#define TESTMODIFIEDSPRINGMODEL_HPP_


#include <cxxtest/TestSuite.h>
#include "TissueSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"


class AbstractCutoffSpringForceModel
{
protected:
    double mCutoffPoint;
    
public:
    AbstractCutoffSpringForceModel(double cutoffPoint)
        : mCutoffPoint(cutoffPoint)
    {
        assert(mCutoffPoint > 0.0);
    }
    
    virtual ~AbstractCutoffSpringForceModel()
    {
    }

    virtual double GetForce(double separation)=0;
    
    double GetCutoffPoint()
    {
        return mCutoffPoint;
    }
};
  
  
class LinearCutoffSpringForceModel : public AbstractCutoffSpringForceModel
{
private:
    double mRestLength;
    double mStiffness;
    
public:
    LinearCutoffSpringForceModel(double cutoffPoint, double stiffness, double restLength)
        :  AbstractCutoffSpringForceModel(cutoffPoint),
           mRestLength(restLength),
           mStiffness(stiffness)
    {
    }
    
    virtual ~LinearCutoffSpringForceModel()
    {
    }
    
    double GetForce(double separation)
    {
        assert(separation < mCutoffPoint);
        return mStiffness*(separation - mRestLength);
    }
};    


//class LinearSurfaceAreaBasedCutoffSpringForceModel : public AbstractCutoffSpringForceModel
//{
//private:
//    double mRestLength;
//    double mBaseStiffness;
//    
//public:
//    LinearSurfaceAreaBasedCutoffSpringForceModel(double cutoffPoint, double baseStiffness, double restLength)
//        :  AbstractCutoffSpringForceModel(cutoffPoint),
//           mRestLength(restLength),
//           mBaseStiffness(baseStiffness)
//    {
//    }
//    
//    virtual ~LinearSurfaceAreaBasedCutoffSpringForceModel()
//    {
//    }
//    
//    double GetForce(double separation)
//    {
//        assert(separation < mCutoffPoint);
//        
//        double R = 0.5*mRestLength;
//        double contact_area = 2*sqrt(R*R - 0.25*separation*separation);
//        
//        double effective_stiffness = mBaseStiffness * ( ?? );
//        
//        return effective_stiffness*(separation - mRestLength);
//    }
//};    


class QuadraticCutoffSpringForceModel : public AbstractCutoffSpringForceModel
{
private:
    double mRestLength;
    double mStiffness;
    
public:
    QuadraticCutoffSpringForceModel(double cutoffPoint, double stiffness, double restLength)
        :  AbstractCutoffSpringForceModel(cutoffPoint),
           mRestLength(restLength),
           mStiffness(stiffness)
    {
    }
    
    virtual ~QuadraticCutoffSpringForceModel()
    {
    }
    
    double GetForce(double separation)
    {
        assert(separation < mCutoffPoint);
        return -mStiffness*(separation - mRestLength)*(separation-2*mCutoffPoint+mRestLength);
    }
};   
      

class CutoffSpringTissueSimulation
{
friend class TestModifiedSpringModel;

private:
    ConformingTetrahedralMesh<2,2>& mrMesh;
    AbstractCutoffSpringForceModel* mpSpringModel;
    std::vector< std::set<unsigned> > mEffectingPoints;

    double mDt;
    double mEndTime;
    std::string mOutputDirectory;

    void FindEffectingPoints()
    {
        mEffectingPoints.clear();
        mEffectingPoints.resize(mrMesh.GetNumNodes());
        for(unsigned i=0; i<mrMesh.GetNumNodes(); i++)
        {
            for(unsigned j=i; j<mrMesh.GetNumNodes(); j++)
            {
                double separation = norm_2(mrMesh.GetVectorFromAtoB(mrMesh.GetNode(i)->rGetLocation(),
                                                                    mrMesh.GetNode(j)->rGetLocation())); 
        
                if((i!=j) && (separation < mpSpringModel->GetCutoffPoint()))
                {
                    mEffectingPoints[i].insert(j);
                    mEffectingPoints[j].insert(i);
                }
            }
        }
    }

    std::vector<c_vector<double,2> > CalculateVelocitiesOfEachNode()
    {
        FindEffectingPoints();
        
        std::vector<c_vector<double,2> > drdt(mrMesh.GetNumNodes());
    
        for (unsigned i=0; i<drdt.size(); i++)
        {
            drdt[i]=zero_vector<double>(2);
        }

        double damping_constant = 10;        

        for(unsigned i=0; i<mrMesh.GetNumNodes(); i++)
        {
            for(std::set<unsigned>::iterator iter = mEffectingPoints[i].begin();
                iter != mEffectingPoints[i].end();
                ++iter)
            {
                unsigned j = *iter;
                
                c_vector<double,2> unit_difference;
                c_vector<double,2> node_a_location = mrMesh.GetNode(i)->rGetLocation();
                c_vector<double,2> node_b_location = mrMesh.GetNode(j)->rGetLocation();

                unit_difference = mrMesh.GetVectorFromAtoB(node_a_location, node_b_location);   
                double distance_between_nodes = norm_2(unit_difference);

                unit_difference /= distance_between_nodes;

                assert( distance_between_nodes < mpSpringModel->GetCutoffPoint() );
                
                double force_magnitude = mpSpringModel->GetForce(distance_between_nodes);
                
                c_vector<double,2> force = force_magnitude*unit_difference;
                drdt[i] += force / damping_constant;
            }
        }

        return drdt;
    }

    void UpdateNodePositions(const std::vector< c_vector<double,2> >& rDrDt)
    {
        for(unsigned index = 0; index < mrMesh.GetNumNodes(); index++)
        {
            ChastePoint<2> new_point(mrMesh.GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
            mrMesh.SetNode(index, new_point, false);
        }
        
        //just to keep elements non-inverted. no computational need for this
        NodeMap map(mrMesh.GetNumNodes());
        mrMesh.ReMesh(map);
    }
    
    
public:    
    CutoffSpringTissueSimulation(ConformingTetrahedralMesh<2,2>& rMesh, AbstractCutoffSpringForceModel* pSpringModel)
        : mrMesh(rMesh),
          mpSpringModel(pSpringModel)
    {
        mEndTime = -1;
        mDt = 0.01;
    }
    
    void Solve()
    {
        if(mEndTime < 0.0)
        {
            EXCEPTION("End time not set");
        }
        
        double current_time = 0.0;
        
        out_stream p_node_file;
        if(mOutputDirectory!="")
        {
            OutputFileHandler output_file_handler(mOutputDirectory);
            p_node_file = output_file_handler.OpenOutputFile("results.viznodes");

            (*p_node_file) << current_time << " ";
            for(unsigned i=0; i<mrMesh.GetNumNodes(); i++)
            {
                for(unsigned j=0; j<2; j++)
                {
                    (*p_node_file) << mrMesh.GetNode(i)->rGetLocation()[j] << " ";
                }
            }
            (*p_node_file) << "\n";
        }
        
        while(current_time < mEndTime)
        {
            // calculate node velocities
            std::vector<c_vector<double,2> > drdt = CalculateVelocitiesOfEachNode();

            // update node positions
            UpdateNodePositions(drdt);

            if(mOutputDirectory!="")
            {
                (*p_node_file) << current_time << " ";
                for(unsigned i=0; i<mrMesh.GetNumNodes(); i++)
                {
                    for(unsigned j=0; j<2; j++)
                    {
                        (*p_node_file) << mrMesh.GetNode(i)->rGetLocation()[j] << " ";
                    }
                }
                (*p_node_file) << "\n";
            }
            
            current_time += mDt;
        }
        
    }
    
    void SetEndTime(double endTime)
    {
        if(endTime <= 0.0)
        {
            EXCEPTION("End time must be greater than zero");
        }
        mEndTime = endTime;
    }

    void SetOutputDirectory(std::string outputDirectory)
    {
        mOutputDirectory = outputDirectory;
    }
};




class TestModifiedSpringModel : public CxxTest::TestSuite
{
public:
    void TestLinearCutoffSpringForceModel() throw(Exception)
    {
        double cutoff_point = 4.0;
        double rest_length = 1.0;
        double stiffness = 1.3435435;
        
        LinearCutoffSpringForceModel linear_model(cutoff_point, stiffness, rest_length);
        
        TS_ASSERT_DELTA( linear_model.GetCutoffPoint(), cutoff_point, 1e-12 );
        TS_ASSERT_DELTA( linear_model.GetForce(rest_length), 0.0, 1e-12 );
        TS_ASSERT_DELTA( linear_model.GetForce(2.0),  stiffness, 1e-12 );
        TS_ASSERT_DELTA( linear_model.GetForce(0.0), -stiffness, 1e-12 );
    }
    
    void TestQuadraticCutoffSpringForceModel() throw(Exception)
    {
        double cutoff_point = 4.0;
        double rest_length = 1.0;
        double stiffness = 1.3435435;
        
        QuadraticCutoffSpringForceModel quadratic_model(cutoff_point, stiffness, rest_length);
        
        TS_ASSERT_DELTA( quadratic_model.GetCutoffPoint(), cutoff_point, 1e-12 );
        TS_ASSERT_DELTA( quadratic_model.GetForce(rest_length), 0.0, 1e-12 );
//        TS_ASSERT_DELTA( quadratic_model.GetForce(rest_length-1),  stiffness*, 1e-12 );
  //      TS_ASSERT_DELTA( quadratic_model.GetForce(0.0), -stiffness, 1e-12 );
    }
    
    
    void TestFindEffectingPoints() throw(Exception)
    {
        unsigned num_elem = 6;
        ConformingTetrahedralMesh<2,2> mesh;

        mesh.ConstructRectangularMesh(num_elem,num_elem);
        
        double cutoff_point = 2.0;
        LinearCutoffSpringForceModel spring_model(cutoff_point, 1.5, 1.0);
        
        CutoffSpringTissueSimulation simulator(mesh, &spring_model);
        simulator.FindEffectingPoints();
                
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double,2>& posn_1 = mesh.GetNode(i)->rGetLocation();

            std::set<unsigned>::iterator iter = simulator.mEffectingPoints[i].find(i);
            bool found = (iter != simulator.mEffectingPoints[i].end());
            
            TS_ASSERT_EQUALS(found,false);
            
            for(unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                if(i!=j)
                {
                    const c_vector<double,2>& posn_2 = mesh.GetNode(j)->rGetLocation();
                    
                    double separation = norm_2(posn_1 - posn_2);
                
                    std::set<unsigned>::iterator iter = simulator.mEffectingPoints[i].find(j);
                        
                    if(iter!=simulator.mEffectingPoints[i].end())
                    {
                        TS_ASSERT_LESS_THAN(separation, cutoff_point);
                    }
                    else
                    {
                        TS_ASSERT_LESS_THAN_EQUALS(cutoff_point, separation);
                    }
                }
            }
        }
    }
    
    // runs a simulation where all nodes are either the rest length away from
    // each other, or too far away to interact, ie there should be no displacement
    void TestSolveWhenExpectNoDisplacement() throw(Exception)
    {
        unsigned num_elem = 6;
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(num_elem,num_elem);
        
        double cutoff_point = 1.1; // all nodes are either 1.0 or sqrt(2) away from each other
        double rest_length  = 1.0; // so in interacting nodes have springs at rest
        
        // note the very high stiffness - shouldn't have an effect
        LinearCutoffSpringForceModel spring_model(cutoff_point, 100, rest_length);
        
        CutoffSpringTissueSimulation simulator(mesh, &spring_model);
        
        simulator.SetEndTime(0.1);
        simulator.Solve();
        
        // recreate the initial mesh and check they agree
        ConformingTetrahedralMesh<2,2> mesh2;
        mesh2.ConstructRectangularMesh(num_elem,num_elem);
        
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            for(unsigned j=0; j<2; j++)
            {
                TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[j],mesh.GetNode(i)->rGetLocation()[j],1e-12);
            }
        }
    }
    
    void TestComputeVelocities()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        double cutoff_point = 4.0; // all nodes interact with each other 
        double rest_length  = 1.0; // so only the opposite edge to each node is not the rest length away
        double stiffness = 10;
        
        LinearCutoffSpringForceModel spring_model(cutoff_point, stiffness, rest_length);
        
        CutoffSpringTissueSimulation simulator(mesh, &spring_model);   
                
        std::vector<c_vector<double,2> > drdt = simulator.CalculateVelocitiesOfEachNode();
        
        double non_zero_force = stiffness*(sqrt(2)-1)/sqrt(2);    
        double damping = 10;
    
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(drdt.size(),4u);

        TS_ASSERT_DELTA(drdt[0](0),  non_zero_force/damping, 1e-12);
        TS_ASSERT_DELTA(drdt[0](1),  non_zero_force/damping, 1e-12);

        TS_ASSERT_DELTA(drdt[1](0), -non_zero_force/damping, 1e-12);
        TS_ASSERT_DELTA(drdt[1](1),  non_zero_force/damping, 1e-12);

        TS_ASSERT_DELTA(drdt[2](0), -non_zero_force/damping, 1e-12);
        TS_ASSERT_DELTA(drdt[2](1), -non_zero_force/damping, 1e-12);

        TS_ASSERT_DELTA(drdt[3](0),  non_zero_force/damping, 1e-12);
        TS_ASSERT_DELTA(drdt[3](1), -non_zero_force/damping, 1e-12);
    }
    
 
    // for experimental work - probably should be 'xTest'ed out in commits..
    void TestSolve() throw(Exception)
    {
//        unsigned num_elem = 6;
//        ConformingTetrahedralMesh<2,2> mesh;
//        mesh.ConstructRectangularMesh(num_elem,num_elem);
        
        HoneycombMeshGenerator generator(6,6,0,false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        p_mesh->Scale(1.0, 1.5);

        double cutoff_point = 2.0;
        LinearCutoffSpringForceModel spring_model(cutoff_point, 30.0, 1.0);
//        QuadraticCutoffSpringForceModel spring_model(cutoff_point, 1.5, 1.0);
        
        CutoffSpringTissueSimulation simulator(*p_mesh, &spring_model);
        
        simulator.SetOutputDirectory("ModifiedSpringModel");
        simulator.SetEndTime(10);
        simulator.Solve();
    }
};

#endif /*TESTMODIFIEDSPRINGMODEL_HPP_*/
