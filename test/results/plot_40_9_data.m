%
% Code to plot the output of the Chaste test
% TestMeinekeLabellingExperiments.hpp
%
% Requires files of the number of cells in each position that were
% labelled in a column, at time = 40 minutes and time = 9 hours.
% 
% Author: Gary Mirams
% Date: 1/8/08

function plot_40_9_data(title_string, this_path)

% This is the number of separate runs of the program
runs = 5;
% This is the number of simulations each run did.
simulations = 50;
% For top of averaged results
crypt_height = 23;

% Don't mess below 'ere!
sections_taken = runs*simulations*2;

t40_results = zeros(100,1);
t9_results = zeros(100,1);
Error40 = zeros(100,1);
Error9 = zeros(100,1);

results_matrix_40 = zeros(100,runs);
results_matrix_9 = zeros(100,runs);

for i=1:runs

    sim_str = num2str(i);
    simulations_str = num2str(simulations);
    
    results_matrix_40(:,i) =  load([this_path '/labelled_t40_totals', simulations_str, '_run_', sim_str, '.dat']);
    results_matrix_9(:,i) = load([this_path '/labelled_t9_totals', simulations_str, '_run_', sim_str, '.dat']);

    t40_results = t40_results + load([this_path '/labelled_t40_totals', simulations_str, '_run_', sim_str, '.dat']);
    t9_results = t9_results + load([this_path '/labelled_t9_totals', simulations_str, '_run_', sim_str, '.dat']);

end

for i=1:100
    Error40(i) = std(results_matrix_40(i,:));
end
for i=1:100
    Error9(i) = std(results_matrix_9(i,:));
end

% Scale to percentages
t40_results = 100.*t40_results./sections_taken;
t9_results = 100.*t9_results./sections_taken;

figure
plot(t9_results,'kx')
hold on
plot(t40_results,'k.')

% Create xlabel
xlabel('Cell Position','FontSize',16);

% Create ylabel
ylabel('% Labelled Cells','FontSize',16);
title(title_string,'FontSize',16);

legend('t = 9:00','t = 0:40')
xlim([0 40])
ylim([0 100])
set(gca,'Fontsize',14)

addpath('../../../../../anim/matlab/'); % Adds the LoadNonConstantLengthData function.

y_all_40min = [];
y_all_9hrs = [];

buckets = 0:1:ceil(crypt_height);

total_num_in_each_bucket_40min = 0*buckets(1:end-1);
total_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

marked_num_in_each_bucket_40min = 0*buckets(1:end-1);
marked_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

vis_nodes_40min = LoadNonConstantLengthData([this_path '/first_lines.txt']);
vis_nodes_9hrs = LoadNonConstantLengthData([this_path '/last_lines.txt']);

if(length(vis_nodes_40min) ~= length(vis_nodes_9hrs))
   error('First and Last line files are of different length');
end;

average_40 = zeros(length(vis_nodes_40min),length(buckets)-1);
average_9 = zeros(length(vis_nodes_40min),length(buckets)-1);

for i=1:length(vis_nodes_40min)

    %	Loop over nodes and if it is marked get y value and plonk in a bucket
    total_num_in_each_bucket_this_exp_40min = 0*buckets(1:end-1);
    total_num_in_each_bucket_this_exp_9hrs = 0*buckets(1:end-1);

    marked_num_in_each_bucket_this_exp_40min = 0*buckets(1:end-1);
    marked_num_in_each_bucket_this_exp_9hrs = 0*buckets(1:end-1);      
    
    num_nodes_40min = (length(vis_nodes_40min{i})-1)/3;
    for j = 1:num_nodes_40min
        y_val = vis_nodes_40min{i}(3*j);
        cell_type = vis_nodes_40min{i}(3*j + 1);
        for k = 1:(length(buckets)-1)
            if y_val >= buckets(k) && y_val < buckets(k+1)
                total_num_in_each_bucket_40min(k) = total_num_in_each_bucket_40min(k) + 1;
                total_num_in_each_bucket_this_exp_40min(k) = total_num_in_each_bucket_this_exp_40min(k) + 1;
                if  cell_type == 5
                    marked_num_in_each_bucket_40min(k) = marked_num_in_each_bucket_40min(k) + 1;
                    marked_num_in_each_bucket_this_exp_40min(k) = marked_num_in_each_bucket_this_exp_40min(k) + 1;
                end

                break;
            end
        end
    end

    
    num_nodes_9hrs = (length(vis_nodes_9hrs{i})-1)/3;
    for j = 1:num_nodes_9hrs
        y_val = vis_nodes_9hrs{i}(3*j);
        cell_type = vis_nodes_9hrs{i}(3*j + 1);
        for k = 1:(length(buckets)-1)
            if y_val >= buckets(k) && y_val < buckets(k+1)
                total_num_in_each_bucket_9hrs(k) = total_num_in_each_bucket_9hrs(k) + 1;
                total_num_in_each_bucket_this_exp_9hrs(k) = total_num_in_each_bucket_this_exp_9hrs(k) + 1;
                if  cell_type == 5
                    marked_num_in_each_bucket_9hrs(k) = marked_num_in_each_bucket_9hrs(k) + 1;
                    marked_num_in_each_bucket_this_exp_9hrs(k) = marked_num_in_each_bucket_this_exp_9hrs(k) + 1;
                end

                break;
            end
        end
    end
    if total_num_in_each_bucket_this_exp_40min ~= 0
        average_40(i,:) = 100*marked_num_in_each_bucket_this_exp_40min./total_num_in_each_bucket_this_exp_40min;
    end
    if total_num_in_each_bucket_this_exp_9hrs ~= 0
        average_9(i,:) = 100*marked_num_in_each_bucket_this_exp_9hrs./total_num_in_each_bucket_this_exp_9hrs;
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

for i=1:length(buckets)-1
    av_error_each_bucket_40(i) = std(average_40(:,i));
    av_error_each_bucket_9(i) = std(average_9(:,i));
end

figure
bar(buckets(1:end-1)+0.5*(buckets(2) - buckets(1)),percent_in_each_bucket_9hrs,'FaceColor',[0.5 0.5 0.5])
hold on
bar(buckets(1:end-1)+0.5*(buckets(2) - buckets(1)),percent_in_each_bucket_40min,'FaceColor',[0.8 0.8 0.8])
xlabel('Height up crypt (cells)','FontSize',16);
ylabel('% Labelled cells','FontSize',16);
title(title_string,'FontSize',16);
ylim([0 100])
xlim([0 20])
legend('t = 9:00','t = 0:40')
set(gca,'Fontsize',14)
hold on
errorbar(buckets(1:end-1)+0.5*(buckets(2) - buckets(1)),percent_in_each_bucket_9hrs,av_error_each_bucket_9,'k.','MarkerSize',1)
hold on
errorbar(buckets(1:end-1)+0.5*(buckets(2) - buckets(1)),percent_in_each_bucket_40min,av_error_each_bucket_40,'k.','MarkerSize',1)

