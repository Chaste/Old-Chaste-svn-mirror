#ifndef DISCRETESYSTEMFORCECALCULATOR_HPP_
#define DISCRETESYSTEMFORCECALCULATOR_HPP_

#include "Meineke2001SpringSystem.hpp"
#include "OutputFileHandler.hpp"


class DiscreteSystemForceCalculator
{
    friend class TestDiscreteSystemForceCalculator;
    
private:
    
    /**
     * Spring system, used to get mesh and calculate forces between nodes.
     */ 
    Meineke2001SpringSystem<2>& mrMeinekeSpringSystem;
    
    /**
     * Small parameter, used in GetSamplingAngles().
     */ 
    double mEpsilon; 
    
    /** The file that the results of CalculateExtremalNormalForces */ 
    out_stream mpStressResultsFile;    
    
    /**
     * Given a node index, returns the set of neighbouring node indices.
     */ 
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        Node<2>* p_node = r_mesh.GetNode(index);
        
        std::set<unsigned> neighbouring_node_indices;
        
        for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element<2,2>* p_element = r_mesh.GetElement(*it);
            for(unsigned i=0; i<p_element->GetNumNodes(); i++)
            {
                unsigned node_index = p_element->GetNodeGlobalIndex(i);
                if(node_index!=index)
                {
                    neighbouring_node_indices.insert(node_index);
                }
            }
        }
        return neighbouring_node_indices;
    }
        
    /**
     * Given a node index and angle of intersecting line in the range (-pi,pi], 
     * returns the tangential and normal forces.
     */ 
    std::vector<double> CalculateFtAndFn(unsigned index, double theta)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
                
        std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(index);
        
        double tangential_force = 0.0;
        double normal_force = 0.0;   
        double alpha;
        
        c_vector<double,2> unit_vec_between_nodes(2);        
        
        for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
            iter != neighbouring_node_indices.end();
            iter++)
        {
            // The method GetAngleBetweenNodes() returns an angle in the range (-pi,pi]
            alpha = r_mesh.GetAngleBetweenNodes(index, *iter);

            assert(alpha <= M_PI);
            assert(alpha > -M_PI);
            
            if ( sin(alpha-theta) > DBL_EPSILON )
            {                
                c_vector<double,2> force_between_nodes = mrMeinekeSpringSystem.CalculateForceBetweenNodes(index,*iter);
                
                unit_vec_between_nodes[0] = cos(alpha);
                unit_vec_between_nodes[1] = sin(alpha);
                
                double plusminus_norm_force = inner_prod(force_between_nodes,unit_vec_between_nodes);         
                tangential_force += plusminus_norm_force * cos(alpha-theta);
                normal_force += plusminus_norm_force * sin(alpha-theta);
            }
        }

        std::vector<double> ret(2);
        ret[0] = tangential_force;
        ret[1] = normal_force;
        
        return ret;
    }
    
    
    /**
     * Given a node index, returns a vector of sampling angles in the range (-pi,pi]
     * that can be used by GetExtremalAngles() to find the locations of local extrema 
     * of the normal force.
     */  
    std::vector<double> GetSamplingAngles(unsigned index)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(index);
        
        std::vector<double> sampling_angles(4*neighbouring_node_indices.size());
        
        unsigned i=0;
        
        for(std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
            iter != neighbouring_node_indices.end();
            iter++)
        {
            // The method GetAngleBetweenNodes() returns an angle in the range (-pi,pi]
            double alpha = r_mesh.GetAngleBetweenNodes(index, *iter);
            
            double alpha_minus_epsilon = alpha - mEpsilon;
            double alpha_plus_epsilon = alpha + mEpsilon;            
            double alpha_plus_pi_minus_epsilon = alpha + M_PI - mEpsilon;
            double alpha_plus_pi_plus_epsilon = alpha + M_PI + mEpsilon;
            
            // Calculate sampling angles in the range (-pi,pi]
            
            #define COVERAGE_IGNORE                    
            if (alpha_minus_epsilon <= -M_PI)
            {
                alpha_minus_epsilon += 2*M_PI;
            }
            #undef COVERAGE_IGNORE
            sampling_angles[i] = alpha_minus_epsilon;  
            
            assert(sampling_angles[i] <= M_PI);
            assert(sampling_angles[i] > -M_PI);
            i++;
            
            if (alpha_plus_epsilon > M_PI)
            {
                alpha_plus_epsilon -= 2*M_PI;
            }
            sampling_angles[i] = alpha_plus_epsilon;
           
            assert(sampling_angles[i] <= M_PI);
            assert(sampling_angles[i] > -M_PI);
            i++;
            
            if (alpha_plus_pi_minus_epsilon > M_PI)
            {
                alpha_plus_pi_minus_epsilon -= 2*M_PI;
            }   
            sampling_angles[i] = alpha_plus_pi_minus_epsilon;
            
            assert(sampling_angles[i] <= M_PI);
            assert(sampling_angles[i] > -M_PI);
            i++;
            
            if (alpha_plus_pi_plus_epsilon > M_PI)
            {
                alpha_plus_pi_plus_epsilon -= 2*M_PI;
            }   
            sampling_angles[i] = alpha_plus_pi_plus_epsilon;
            
            assert(sampling_angles[i] <= M_PI);
            assert(sampling_angles[i] > -M_PI);
            i++;            
        }
        
        sort(sampling_angles.begin(), sampling_angles.end());        
        return sampling_angles;
    }
    
    
    /**
     * Given a node index and two sampling angles, finds the location of 
     * the root of the tangential force in the interval between the two 
     * angles. There is no guarantee that this will lie in (-pi,pi].
     */ 
    double GetLocalExtremum(unsigned index, double angle1, double angle2)
    {   
        // We always pass in angle1 and angle2 such that angle1<angle2, 
        // but note that angle1 may be <M_PI           
        assert(angle1 < angle2);
        
        double tolerance = 1e-5;
        unsigned counter = 0;
                
        double previous_angle;        
        double current_error;     
        double current_angle = angle1; 
        
        current_error = angle2 - angle1;            
        std::vector<double> current_ft_and_fn(2);
        
        while (current_error > tolerance)
        {
            previous_angle = current_angle;
            current_ft_and_fn = CalculateFtAndFn(index, current_angle);
            current_angle -= current_ft_and_fn[0]/current_ft_and_fn[1];
            current_error = fabs(current_angle - previous_angle);
            counter++;
        }
        
        assert( current_angle>angle1 && current_angle<angle2 );            
        assert( current_error < tolerance );        
        
        return current_angle;        
    }
    

    /**
     * Given a vector of sampling angles in the range (-pi,pi], returns a vector 
     * of extremal angles, i.e. angles at which local extrema of the normal force 
     * occur, again in the range (-pi,pi].
     */ 
    std::vector<double> GetExtremalAngles(unsigned index, std::vector<double> samplingAngles)
    {
        std::vector<double> extremal_angles;
        std::vector<double> ft_and_fn(2);
        std::vector<double> tangential_force(samplingAngles.size());
        
        for (unsigned i=0; i<samplingAngles.size(); i++)
        {
            ft_and_fn = CalculateFtAndFn(index, samplingAngles[i]);
            tangential_force[i] = ft_and_fn[0];
        }         
        
        unsigned n = samplingAngles.size()-1;
        
        for (unsigned i=0; i<n; i++)
        {
            if ( ( tangential_force[i%n]>0 && tangential_force[(i+1)%n]<0 ) ||  
                 ( tangential_force[i%n]<0 && tangential_force[(i+1)%n]>0 ) )
            {
                double next_extremal_angle;
                
                // If we are in the interval that crosses the branch line at pi, 
                // then subtract 2*pi from the positive angle 
                if (i==n-1)
                {
                    samplingAngles[i%n] -= 2*M_PI;
                }
                
                if (samplingAngles[(i+1)%n] - samplingAngles[i%n] < 2*mEpsilon + 1e-6 )
                {
                    // If we find a jump through zero, then the local extremum is
                    // simply at the mid-point of the interval
                    next_extremal_angle = (samplingAngles[(i+1)%n] + samplingAngles[i%n])/2.0;
                }
                else
                {
                    // Otherwise we need to find it using Newton's method
                    next_extremal_angle = GetLocalExtremum(index, samplingAngles[i%n], samplingAngles[(i+1)%n]);    
                }
                
                if (next_extremal_angle <= -M_PI)
                {
                    next_extremal_angle += 2*M_PI;
                }
                assert(next_extremal_angle>-M_PI && next_extremal_angle<=M_PI); 
                extremal_angles.push_back(next_extremal_angle);
            }
        }
        
        return extremal_angles;        
    }
    
public:
    
    DiscreteSystemForceCalculator(Meineke2001SpringSystem<2>& rMeinekeSpringSystem)
        : mrMeinekeSpringSystem(rMeinekeSpringSystem),
          mEpsilon(0.01)
    {
    }
    
    std::vector< std::vector<double> > CalculateExtremalNormalForces()
    {
        std::vector< std::vector<double> > extremal_normal_forces;
        
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        
        std::vector<double> minimum_normal_forces(r_mesh.GetNumNodes());
        std::vector<double> maximum_normal_forces(r_mesh.GetNumNodes());
                
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            std::vector<double> sampling_angles = GetSamplingAngles(i);            
            std::vector<double> extremal_angles = GetExtremalAngles(i,sampling_angles);

            double minimum_normal_force_for_node_i = DBL_MAX;
            double maximum_normal_force_for_node_i = -DBL_MAX;
            
            double current_normal_force;
            
            for (unsigned j=0; j<extremal_angles.size(); j++)
            {
                current_normal_force = CalculateFtAndFn(i,extremal_angles[j])[1];
                if (current_normal_force > maximum_normal_force_for_node_i)
                {
                    maximum_normal_force_for_node_i = current_normal_force;
                }
                
                if (current_normal_force < minimum_normal_force_for_node_i)
                {
                    minimum_normal_force_for_node_i = current_normal_force;
                }
            }   
            
            minimum_normal_forces[i] = minimum_normal_force_for_node_i;
            maximum_normal_forces[i] = maximum_normal_force_for_node_i;            
        } 
                
        extremal_normal_forces.push_back(minimum_normal_forces);
        extremal_normal_forces.push_back(maximum_normal_forces);
        
        return extremal_normal_forces;
    }
    
    void WriteResultsToFile(std::string simulationOutputDirectory)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetDimensionalisedTime();
        std::ostringstream time_string;
        time_string << time;
        std::string results_directory = simulationOutputDirectory +"/results_from_time_" + time_string.str();
            
        OutputFileHandler output_file_handler2(results_directory+"/vis_results/",false); 
        mpStressResultsFile = output_file_handler2.OpenOutputFile("results.vizstress");
        
        (*mpStressResultsFile) <<  time << "\t";
        
        double global_index;
        double x;
        double y;
        double minimum;
        double maximum;

        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        
        std::vector< std::vector<double> > extremal_normal_forces = CalculateExtremalNormalForces();
               
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            global_index = (double) i;
            
            x = r_mesh.GetNode(i)->rGetLocation()[0];
            y = r_mesh.GetNode(i)->rGetLocation()[1];
            
            minimum = extremal_normal_forces[0][i]; 
            maximum = extremal_normal_forces[1][i];
            (*mpStressResultsFile) << global_index << " " << x << " " << y << " " << minimum << " " << maximum << " ";
        }
        
        (*mpStressResultsFile) << "\n";        
        mpStressResultsFile->close();
    }

};
#endif /*DISCRETESYSTEMFORCECALCULATOR_HPP_*/
