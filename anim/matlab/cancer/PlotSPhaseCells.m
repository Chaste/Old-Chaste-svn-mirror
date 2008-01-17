% PlotCellsInSPhase

Start_time = 300;
Num_experiments = 50;
crypt_height = 30;

y_all_40min = [];
y_all_9hrs = [];

buckets = 0:1:ceil(crypt_height);
total_num_in_each_bucket_40min = 0*buckets(1:end-1);
total_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

marked_num_in_each_bucket_40min = 0*buckets(1:end-1);
marked_num_in_each_bucket_9hrs = 0*buckets(1:end-1);

%	Get data for each exp
for i=1:Num_experiments
	disp('')
	Experiment_time = Start_time + 10*(i-1);
	FileName = ['/local/pmxaw/MeinekeLabellingExperiment/results_from_time_' int2str(Experiment_time) '.667/vis_results/first_line.txt'];
	
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
	
	
	FileName = ['/local/pmxaw/MeinekeLabellingExperiment/results_from_time_' int2str(Experiment_time) '.667/vis_results/last_line.txt'];
	
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
subplot(2,1,1)
bar(buckets(1:end-1)+0.5,percent_in_each_bucket_40min)
subplot(2,1,2)
bar(buckets(1:end-1)+0.5,percent_in_each_bucket_9hrs)
