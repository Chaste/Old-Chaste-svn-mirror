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
title_string = 'Simple Wnt Cells in Sunter iii) Geometry';
crypt_height = 30;
path = '/local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter3/';
% End of setup

addpath('../');	% Adds the LoadNonConstantLengthData function.

y_all_40min = [];
y_all_9hrs = [];

buckets = 0:1:ceil(crypt_height);

total_num_in_each_bucket_40min = 0*buckets(1:end-1);
total_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

marked_num_in_each_bucket_40min = 0*buckets(1:end-1);
marked_num_in_each_bucket_9hrs = 0*buckets(1:end-1);


vis_nodes_40min = LoadNonConstantLengthData([path 'first_lines.txt']);
vis_nodes_9hrs = LoadNonConstantLengthData([path 'last_lines.txt']);

if(length(vis_nodes_40min) ~= length(vis_nodes_9hrs))
   error('First and Last line files are of different length');
end;

for i=1:length(vis_nodes_40min)

    %	Loop over nodes and if it is marked get y value and plonk in a bucket
    
    num_nodes_40min = (length(vis_nodes_40min{i})-1)/3;
    for j = 1:num_nodes_40min
        y_val = vis_nodes_40min{i}(3*j);
        cell_type = vis_nodes_40min{i}(3*j + 1);
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

    num_nodes_9hrs = (length(vis_nodes_9hrs{i})-1)/3;
    for j = 1:num_nodes_9hrs
        y_val = vis_nodes_9hrs{i}(3*j);
        cell_type = vis_nodes_9hrs{i}(3*j + 1);
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
bar(buckets(1:end-1)+0.5*(buckets(2) - buckets(1)),percent_in_each_bucket_9hrs,'r')
hold on
bar(buckets(1:end-1)+0.5*(buckets(2) - buckets(1)),percent_in_each_bucket_40min,'b')
title([title_string '. After 40 minutes and 9 hours for ' num2str(length(vis_nodes_40min)) ' experiments.']);
xlabel('Height up crypt (cells)');
ylabel('% of labelled cells');
ylim([0 100]);
