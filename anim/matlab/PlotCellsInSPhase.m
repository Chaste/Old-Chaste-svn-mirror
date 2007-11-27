% PlotCellsInSPhase.m
%
% This code needs FirstAndLastScript.m to have been run on the relevant set
% of results.
%
% It reads in the positions of all cells at the beginning and end of a 
% Meineke-stlye labelling experiment and plots the percentages of cells 
% that are labelled in ranges 
%

close all
clear

% Experiment Setup
title_string = 'Meineke Cells in Sunter i) Geometry';
runs = 5;
Start_time = 300;
Num_experiments = 51;
crypt_height = 30;
file_path(1,:) = '/local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter1/2007-11-24-18-10/MeinekeLabellingExperiment';
file_path(2,:) = '/local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter1/2007-11-24-18-12/MeinekeLabellingExperiment';
file_path(3,:) = '/local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter1/2007-11-24-18-15/MeinekeLabellingExperiment';
file_path(4,:) = '/local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter1/2007-11-24-18-16/MeinekeLabellingExperiment';
file_path(5,:) = '/local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter1/2007-11-24-18-19/MeinekeLabellingExperiment';
% End of setup

y_all_40min = [];
y_all_9hrs = [];

buckets = 0:1:ceil(crypt_height);
total_num_in_each_bucket_40min = 0*buckets(1:end-1);
total_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

marked_num_in_each_bucket_40min = 0*buckets(1:end-1);
marked_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

for run = 1:runs
    %	Get data for each exp
    for i=1:Num_experiments
        disp('')
        Experiment_time = Start_time + 10*(i-1);
        temp_string = deblank(file_path(run,:));
        FileName = [temp_string '/results_from_time_' int2str(Experiment_time) '.667/vis_results/first_line.txt'];

        vis_nodes_40min = LoadNonConstantLengthData(FileName);

        %	Loop over nodes and if it is marked get y value and plonk in a bucket
        num_nodes_40min = (length(vis_nodes_40min{1})-1)/3;

        for j = 1:num_nodes_40min
            y_val = vis_nodes_40min{1}(3*j);
            cell_type = vis_nodes_40min{1}(3*j + 1);
            for k = 1:length(buckets)
                if y_val >= buckets(k) && y_val < buckets(k+1)
                    total_num_in_each_bucket_40min(k) = total_num_in_each_bucket_40min(k) + 1;
                    if  cell_type == 5
                        marked_num_in_each_bucket_40min(k) = marked_num_in_each_bucket_40min(k) + 1;
                    end

                    break;
                end
            end
        end


        FileName = [temp_string '/results_from_time_' int2str(Experiment_time) '.667/vis_results/last_line.txt'];

        vis_nodes_9hrs = LoadNonConstantLengthData(FileName);

        num_nodes_9hrs = (length(vis_nodes_9hrs{1})-1)/3;


        for j = 1:num_nodes_9hrs
            y_val = vis_nodes_9hrs{1}(3*j);
            cell_type = vis_nodes_9hrs{1}(3*j + 1);
            for k = 1:length(buckets)
                if y_val >= buckets(k) && y_val < buckets(k+1)
                    total_num_in_each_bucket_9hrs(k) = total_num_in_each_bucket_9hrs(k) + 1;
                    if  cell_type == 5
                        marked_num_in_each_bucket_9hrs(k) = marked_num_in_each_bucket_9hrs(k) + 1;
                    end

                    break;
                end
            end
        end

    end

end

percent_in_each_bucket_40min = 0*total_num_in_each_bucket_40min;
percent_in_each_bucket_9hrs = 0*total_num_in_each_bucket_9hrs;

for i = 1:length(buckets)-1
	if total_num_in_each_bucket_40min(i) ~= 0
		percent_in_each_bucket_40min(i) = 100*marked_num_in_each_bucket_40min(i) / total_num_in_each_bucket_40min(i);
	end
	if total_num_in_each_bucket_9hrs(i) ~= 0
		percent_in_each_bucket_9hrs(i) = 100*marked_num_in_each_bucket_9hrs(i) / total_num_in_each_bucket_9hrs(i);
	end
end

figure;
bar(buckets(1:end-1)+0.5,percent_in_each_bucket_9hrs,'r')
hold on
bar(buckets(1:end-1)+0.5,percent_in_each_bucket_40min,'b')
title([title_string '. After 40 minutes and 9 hours for ' int2str(runs*Num_experiments) ' experiments.']);
xlabel('Height up crypt (cells)');
ylabel('% of labelled cells');
ylim([0 100]);
