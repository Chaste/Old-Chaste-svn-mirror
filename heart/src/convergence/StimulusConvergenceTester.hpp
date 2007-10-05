#ifndef STIMULUSCONVERGENCETESTER_HPP_
#define STIMULUSCONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class StimulusConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    double ExpectedStimulus;
    double Increment;
    unsigned FirstMesh;
    unsigned PlotPoints;
    void PopulateStandardResult(double standardResult[])
    {
        assert(this->PopulatedResult==false);
  
        this->MeshNum=FirstMesh;
        if (DIM==2 || DIM==3)
        {
	      HardCodedResult(standardResult);
          ExpectedStimulus = this->AbsoluteStimulus;
        
        } 
	else
	{
        CuboidMeshConstructor<DIM> constructor;
        std::string mesh_pathname;
 
        mesh_pathname = constructor.Construct(this->MeshNum);
        
        
        assert (this->UseAbsoluteStimulus==true);
        GeneralPlaneStimulusCellFactory<CELL, DIM> cell_factory(this->OdeTimeStep, this->AbsoluteStimulus, true);
        CARDIAC_PROBLEM cardiac_problem(&cell_factory);
        
        cardiac_problem.SetMeshFilename(mesh_pathname);
        cardiac_problem.SetOutputDirectory ("Convergence");
        cardiac_problem.SetOutputFilenamePrefix ("Results");
        
        cardiac_problem.SetEndTime(simulation_time);   // ms
        cardiac_problem.SetLinearSolverRelativeTolerance(this->KspRtol);

        cardiac_problem.SetPdeTimeStep(this->PdeTimeStep);
        
        assert(fabs(0.04/this->PdeTimeStep - round(0.04/this->PdeTimeStep)) <1e-15 );
        cardiac_problem.SetPrintingTimeStep(0.04);  //Otherwise we can't take the timestep down to machine precision without generating thousands of output files
        cardiac_problem.Initialise();
        
        this->DisplayRun();
        try
        {
            cardiac_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout<<"Warning - this run threw an exception.  Check convergence results\n"; 
            assert(0); //Can't continue                
        }
        // Calculate positions of nodes 1/4 and 3/4 through the mesh
        unsigned third_quadrant_node;
        unsigned first_quadrant_node;
        switch(DIM)
        {
            case 1:
            {
                first_quadrant_node = (unsigned) (0.25*constructor.NumElements);
                third_quadrant_node = (unsigned) (0.75*constructor.NumElements);
                assert(cardiac_problem.rGetMesh().GetNode(first_quadrant_node)->rGetLocation()[0]==0.25*mesh_width);
                assert(cardiac_problem.rGetMesh().GetNode(third_quadrant_node)->rGetLocation()[0]==0.75*mesh_width);
                break;
            }
            case 2:
            {
                unsigned n= (unsigned) pow (2, this->MeshNum+2);
                first_quadrant_node =   (n+1)*(n/2)+  n/4 ;
                third_quadrant_node =   (n+1)*(n/2)+3*n/4 ;
                break;
            }
            case 3:
            {
                const unsigned first_quadrant_nodes_3d[5]={61, 362, 2452, 17960, 137296};
                const unsigned third_quadrant_nodes_3d[5]={63, 366, 2460, 17976, 137328};
                assert(this->MeshNum<5);
                first_quadrant_node = first_quadrant_nodes_3d[this->MeshNum];
                third_quadrant_node = third_quadrant_nodes_3d[this->MeshNum];
                break;
            }
            
            default:
                assert(0);
        }
        
        Node<DIM>* fqn = cardiac_problem.rGetMesh().GetNode(first_quadrant_node);
        Node<DIM>* tqn = cardiac_problem.rGetMesh().GetNode(third_quadrant_node);
        assert(fqn->rGetLocation()[0]==0.25*mesh_width);
        assert(tqn->rGetLocation()[0]==0.75*mesh_width);
        for (unsigned coord=1; coord<DIM; coord++)
        {
            assert(fqn->rGetLocation()[coord]==0.5*mesh_width);
            assert(tqn->rGetLocation()[coord]==0.5*mesh_width);
        }
        
        OutputFileHandler results_handler("Convergence", false);
        ColumnDataReader results_reader(results_handler.GetTestOutputDirectory(), "Results", false);
        
        
        // Write out the time series for the node at first and third quadrant
        {
            std::vector<double> transmembrane_potential=results_reader.GetValues("V", third_quadrant_node);
            std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
            OutputFileHandler plot_file_handler("ConvergencePlots", false);
            std::stringstream plot_file_name_stream;
            plot_file_name_stream<< "Node1_"<< "base" << "_timestep.csv";
            out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
            for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
            {
                (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
            }
            p_plot_file->close();
            
            
            for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
            {
                standardResult[data_point] = transmembrane_potential[data_point];
            }
        }
        
        {
            std::vector<double> transmembrane_potential=results_reader.GetValues("V", first_quadrant_node);
            std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
            OutputFileHandler plot_file_handler("ConvergencePlots", false);
            std::stringstream plot_file_name_stream;
            plot_file_name_stream<< "Node2_"<< "base" << "_timestep.csv";
            out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
            for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
            {
                (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
            }
            p_plot_file->close();
        }                    
        
        this->MeshNum++;
        ExpectedStimulus = this->AbsoluteStimulus*2;//For next refinement
        
	}
        
        this->PopulatedResult=true;
        Increment=-ExpectedStimulus/8;
        PlotPoints=8;
        this->AbsoluteStimulus = ExpectedStimulus-PlotPoints*Increment/2;
        
       
    }
    
    void SetInitialConvergenceParameters()
    {
        //this->MeshNum=1; //This is set later using FirstMesh
        this->FixedResult=true;
        this->UseAbsoluteStimulus=true;
        
        double mesh_size=pow(2,FirstMesh+2);
        switch (DIM)
        {
            case 1:
            {
                this->AbsoluteStimulus=-1e7*mesh_size/64.0;
                //wrt mesh 4 which has 64 elements in 1D
                break;
            }
            case 2:
            {
                this->AbsoluteStimulus=-1e7 * mesh_size / 64.0;



                break;
            }
            default:
            {
                this->AbsoluteStimulus=-1e7 * mesh_size / 64.0;
            }
        }
        
        //Real picky for now
        this->RelativeConvergenceCriterion=-1; //Never converge
    }
    
    void UpdateConvergenceParameters()
    {
        this->AbsoluteStimulus += Increment;
    }
    
    bool GiveUpConvergence()
    {
        return (this->AbsoluteStimulus>ExpectedStimulus+PlotPoints*Increment/2); //3.5e7
    }
    double Abscissa()
    {
        return this->AbsoluteStimulus;
    }


    void HardCodedResult(double standardResult[])
    {
      /* This is the result of the "converged" test in 1-D
       *  Solving with a space step of 0.0015625 cm (mesh 5)
       *  Solving with a time step of 0.005 ms
       *  Solving with an ode time step of 0.0025 ms
       *  Solving with a KSP relative tolerance of 1e-08
       *  (Stimulus = -2e7)
       */
      std::cout<<"Standard result copied from 1d\n";

        standardResult[0] = -83.853;
        standardResult[1] = -83.853;
        standardResult[2] = -83.853;
        standardResult[3] = -83.853;
        standardResult[4] = -83.853;
        standardResult[5] = -83.854;
        standardResult[6] = -83.854;
        standardResult[7] = -83.854;
        standardResult[8] = -83.854;
        standardResult[9] = -83.854;
        standardResult[10] = -83.854;
        standardResult[11] = -83.854;
        standardResult[12] = -83.854;
        standardResult[13] = -83.854;
        standardResult[14] = -83.854;
        standardResult[15] = -83.853;
        standardResult[16] = -83.852;
        standardResult[17] = -83.85;
        standardResult[18] = -83.847;
        standardResult[19] = -83.842;
        standardResult[20] = -83.835;
        standardResult[21] = -83.825;
        standardResult[22] = -83.811;
        standardResult[23] = -83.793;
        standardResult[24] = -83.769;
        standardResult[25] = -83.739;
        standardResult[26] = -83.7;
        standardResult[27] = -83.65;
        standardResult[28] = -83.589;
        standardResult[29] = -83.513;
        standardResult[30] = -83.419;
        standardResult[31] = -83.304;
        standardResult[32] = -83.163;
        standardResult[33] = -82.993;
        standardResult[34] = -82.787;
        standardResult[35] = -82.537;
        standardResult[36] = -82.236;
        standardResult[37] = -81.873;
        standardResult[38] = -81.435;
        standardResult[39] = -80.908;
        standardResult[40] = -80.274;
        standardResult[41] = -79.511;
        standardResult[42] = -78.595;
        standardResult[43] = -77.493;
        standardResult[44] = -76.17;
        standardResult[45] = -74.581;
        standardResult[46] = -72.674;
        standardResult[47] = -70.388;
        standardResult[48] = -67.648;
        standardResult[49] = -64.366;
        standardResult[50] = -60.44;
        standardResult[51] = -55.748;
        standardResult[52] = -50.149;
        standardResult[53] = -43.497;
        standardResult[54] = -35.686;
        standardResult[55] = -26.783;
        standardResult[56] = -17.232;
        standardResult[57] = -7.9103;
        standardResult[58] = 0.24791;
        standardResult[59] = 6.7448;
        standardResult[60] = 11.602;
        standardResult[61] = 15.113;
        standardResult[62] = 17.611;
        standardResult[63] = 19.379;
        standardResult[64] = 20.625;
        standardResult[65] = 21.502;
        standardResult[66] = 22.118;
        standardResult[67] = 22.549;
        standardResult[68] = 22.854;
        standardResult[69] = 23.073;
        standardResult[70] = 23.24;
        standardResult[71] = 23.38;
        standardResult[72] = 23.518;
        standardResult[73] = 23.67;
        standardResult[74] = 23.846;
        standardResult[75] = 24.044;
        standardResult[76] = 24.253;
        standardResult[77] = 24.459;
        standardResult[78] = 24.649;
        standardResult[79] = 24.815;
        standardResult[80] = 24.952;
        standardResult[81] = 25.057;
        standardResult[82] = 25.13;
        standardResult[83] = 25.173;
        standardResult[84] = 25.187;
        standardResult[85] = 25.174;
        standardResult[86] = 25.137;
        standardResult[87] = 25.078;
        standardResult[88] = 25;
        standardResult[89] = 24.903;
        standardResult[90] = 24.792;
        standardResult[91] = 24.667;
        standardResult[92] = 24.53;
        standardResult[93] = 24.383;
        standardResult[94] = 24.227;
        standardResult[95] = 24.064;
        standardResult[96] = 23.893;
        standardResult[97] = 23.718;
        standardResult[98] = 23.537;
        standardResult[99] = 23.353;
        standardResult[100] = 23.165;
        standardResult[101] = 22.974;
        standardResult[102] = 22.782;
        standardResult[103] = 22.587;
        standardResult[104] = 22.391;
        standardResult[105] = 22.194;
        standardResult[106] = 21.997;
        standardResult[107] = 21.799;
        standardResult[108] = 21.601;
        standardResult[109] = 21.403;
        standardResult[110] = 21.205;
        standardResult[111] = 21.007;
        standardResult[112] = 20.81;
        standardResult[113] = 20.614;
        standardResult[114] = 20.418;
        standardResult[115] = 20.224;
        standardResult[116] = 20.03;
        standardResult[117] = 19.837;
        standardResult[118] = 19.645;
        standardResult[119] = 19.455;
        standardResult[120] = 19.265;
        standardResult[121] = 19.077;
        standardResult[122] = 18.891;
        standardResult[123] = 18.705;
        standardResult[124] = 18.521;
        standardResult[125] = 18.338;
        standardResult[126] = 18.157;
        standardResult[127] = 17.977;
        standardResult[128] = 17.799;
        standardResult[129] = 17.621;
        standardResult[130] = 17.446;
        standardResult[131] = 17.272;
        standardResult[132] = 17.099;
        standardResult[133] = 16.928;
        standardResult[134] = 16.758;
        standardResult[135] = 16.59;
        standardResult[136] = 16.424;
        standardResult[137] = 16.258;
        standardResult[138] = 16.095;
        standardResult[139] = 15.933;
        standardResult[140] = 15.772;
        standardResult[141] = 15.613;
        standardResult[142] = 15.455;
        standardResult[143] = 15.299;
        standardResult[144] = 15.144;
        standardResult[145] = 14.991;
        standardResult[146] = 14.839;
        standardResult[147] = 14.689;
        standardResult[148] = 14.54;
        standardResult[149] = 14.392;
        standardResult[150] = 14.246;
        standardResult[151] = 14.102;
        standardResult[152] = 13.959;
        standardResult[153] = 13.817;
        standardResult[154] = 13.677;
        standardResult[155] = 13.538;
        standardResult[156] = 13.4;
        standardResult[157] = 13.264;
        standardResult[158] = 13.129;
        standardResult[159] = 12.996;
        standardResult[160] = 12.864;
        standardResult[161] = 12.733;
        standardResult[162] = 12.604;
        standardResult[163] = 12.476;
        standardResult[164] = 12.349;
        standardResult[165] = 12.224;
        standardResult[166] = 12.099;
        standardResult[167] = 11.977;
        standardResult[168] = 11.855;
        standardResult[169] = 11.735;
        standardResult[170] = 11.616;
        standardResult[171] = 11.498;
        standardResult[172] = 11.381;
        standardResult[173] = 11.266;
        standardResult[174] = 11.152;
        standardResult[175] = 11.039;
        standardResult[176] = 10.927;
        standardResult[177] = 10.816;
        standardResult[178] = 10.707;
        standardResult[179] = 10.599;
        standardResult[180] = 10.492;
        standardResult[181] = 10.386;
        standardResult[182] = 10.281;
        standardResult[183] = 10.177;
        standardResult[184] = 10.074;
        standardResult[185] = 9.9729;
        standardResult[186] = 9.8725;
        standardResult[187] = 9.7732;
        standardResult[188] = 9.6749;
        standardResult[189] = 9.5777;
        standardResult[190] = 9.4815;
        standardResult[191] = 9.3864;
        standardResult[192] = 9.2923;
        standardResult[193] = 9.1992;
        standardResult[194] = 9.1072;
        standardResult[195] = 9.0161;
        standardResult[196] = 8.926;
        standardResult[197] = 8.8369;
        standardResult[198] = 8.7488;
        standardResult[199] = 8.6617;
        standardResult[200] = 8.5754;
    }

};

#endif /*STIMULUSCONVERGENCETESTER_HPP_*/
