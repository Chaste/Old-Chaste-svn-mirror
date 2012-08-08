function script_for_particular_measures
%SCRIPT Summary of this function goes here
%   Detailed explanation goes here


%% Settings to play with
fit_type = 'linear';
%fit_type = 'quadratic';
whose_code = @my_classification_code; % My manual classification
%whose_code = @classification_code; % Matlab's classify() function

%make_new_variables = true;  % When reducing dimension we recast the problem
                            % into new more discriminant variables.
make_new_variables = false; % When reducing dimension we pick the 
                            % best of the original variables

% Setup
dimension_of_analysis = 1;                            
plotting_on = true;
combine_cats_1_and_2 = true;
remove_ranolazine = true

measure1 = 316; % 8 = hERG IC50, 11 = Redfern measure, 316 = Grandi APD90 Max Effect EFTPC
measure2 = 17;

if dimension_of_analysis==1
    plotting_on = false;
end

%% Start of code - mess below here with caution!
cd ..
data = load_drug_data();
cd classification

% Get the redfern categories to use as the discriminant 'groups'.
output_groups = data.drug_redfern;

%% Load up all possible inputs
load('compiled_results.mat');
original_data = compiled_results.measures;
original_names = compiled_results.measure_names';
for i=1:length(original_names)
   fprintf('%i \t %s \n',i,original_names{i})
end

%% Choose Different drugs with the following inputs available in the dataset
%%% HergIC50 and dose
%use_this_data = find(abs(original_data(:,3) + 1)>1e-6 & original_data(:,5)>0);

%%% All three IC50 values
%use_this_data = find(abs(original_data(:,1) + 1)>1e-6 & abs(original_data(:,2) + 1)>1e-6 & abs(original_data(:,3) + 1)>1e-6);

%% All three IC50 values and dose data
use_this_data = find(abs(original_data(:,1) + 1)>1e-6 & abs(original_data(:,2) + 1)>1e-6 & abs(original_data(:,3) + 1)>1e-6 & original_data(:,4)>0);

if combine_cats_1_and_2
    idxes = find(data.drug_redfern(:)==1 | data.drug_redfern(:)==2);
    output_groups(idxes) = 2;
    data.drug_redfern(idxes)=2;
end

if remove_ranolazine
    use_this_data = use_this_data(use_this_data~=12);
end

% Select only this data
input_data = original_data(use_this_data,:);
redfern_cats = data.drug_redfern(use_this_data);
output_groups = output_groups(use_this_data);

%% Just display the drugs we are going to use with their category...
fprintf('Running classification for %i drugs:\n',length(use_this_data))
for i=1:length(use_this_data)
    fprintf('Drug = %s,\t Redfern Category = %i\n',data.drug_names{use_this_data(i)},output_groups(i))
end
fprintf('\n')

% This snippet is useful for choosing exactly which inputs to plot...

%% CHOOSE WHICH INPUTS TO USE FOR DISCRIMINATION

%% Make training and validation data N dimensional
use_these_vars = 1:dimension_of_analysis;

% % Redfern option modifies this and uses one dimensional discrimination
% input_choice_data = [input_data(:,11) input_data(:,8) ];
% labels = {original_names{11},original_names{8} };
% title_text = 'Redfern Classification';

% % Yi's suggestions are here
% input_choice_data = [input_data(:,9) input_data(:,10)];
% labels = {original_names{9} , original_names{10}};
% title_text = 'GSK Classification';

% Our first suggestion here
input_choice_data = [input_data(:,measure1) input_data(:,measure2)];
labels = {original_names{measure1}, original_names{measure2} };
title_text = 'Simulated Marker Classification';

% This has to be here in case the operations above reorder the
% variables.
input_choice_data = input_choice_data(:,use_these_vars);

%% Run leave-one-out classification analysis
error_score = 0;
for run=1:length(use_this_data)
    [v_score] = feval(whose_code,run, input_choice_data, output_groups,fit_type,plotting_on);
    error_score(run) = v_score;
    figure(1)
    xlabel(labels{1})
    ylabel(labels{2})
end
title(title_text)

name = labels{1};
if dimension_of_analysis==2
    name = [name ' and ' labels{2}];
end

% Make our own histograms so that we can control the size of the bars
% properly... hist() mucks them up...
buckets = zeros(5,1);
for i=1:length(error_score)
    if abs(error_score(i))<1e-6
        buckets(1) = buckets(1)+1;
    elseif (abs(error_score(i)-1)<1e-6)
        buckets(2) = buckets(2)+1;
    elseif  (abs(error_score(i)-2)<1e-6)
        buckets(3) = buckets(3)+1;
    elseif  (abs(error_score(i)-3)<1e-6)
        buckets(4) = buckets(4)+1;
    elseif  (abs(error_score(i)-4)<1e-6)
        buckets(5) = buckets(5)+1;
    end
end
assert(sum(buckets)==length(error_score));

% This score is the mean difference between the predicted value and the
% actual value. It will work for any discrimination technique.
figure
bar([0 1 2 3 4],buckets,'BarWidth',0.5);
%hist(error_score)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k')

xlabel('Error in categorization','FontSize',18)
ylabel('Number of drugs','FontSize',18)
ylim([0 23])
xlim([-0.5 max(redfern_cats)-min(redfern_cats)+0.5])
set(gca,'XTick',0:4,'FontSize',16)
mean_error_in_category = mean(error_score)
std_dev_error_in_category = std(error_score)
name
title(['mean = ' num2str(mean_error_in_category,'% 4.3f')...
    ', std deviation = ' num2str(std_dev_error_in_category,'% 4.3f')],'FontSize',18)

%title([name ', mean = ' num2str(mean_error_in_category)...
%    ', std deviation = ' num2str(std_dev_error_in_category)],'FontSize',18)

% figure(2)
% scatter(-training,-validation,'rx')
% xlabel('Training error (- log likelihood)')
% ylabel('Validation error (- log likelihood)')



