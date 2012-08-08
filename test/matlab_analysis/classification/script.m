function script

close all

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
dimension_of_analysis =2; % When changing this remember to comment in/out an 'if' in the loops below.
plotting_on = false;
remove_cat_1 = false;
remove_ranolazine = true;

combine_cats_1_and_2 = true;

%% Start of code - mess below 'ere with caution!
cd ..
data = load_drug_data();
cd classification

%% Load up all possible inputs
load('compiled_results.mat');
original_data = compiled_results.measures;
original_names = compiled_results.measure_names';
assert(length(original_data)==length(original_names));
% for i=1:length(original_names)
%    fprintf('%i & \t %s \\\\ \n',i,original_names{i})
% end
% pause
fprintf('Total = %i measures.\n',length(original_names));

%% Choose Different drugs with the following inputs available in the dataset
%%% HergIC50 and dose
%use_this_data = find(abs(original_data(:,3) + 1)>1e-6 & original_data(:,5)>0);

%%% All three IC50 values
%use_this_data = find(abs(original_data(:,1) + 1)>1e-6 & abs(original_data(:,2) + 1)>1e-6 & abs(original_data(:,3) + 1)>1e-6);

%% All three IC50 values and dose data
if remove_cat_1
    use_this_data = find(abs(original_data(:,1) + 1)>1e-6 & abs(original_data(:,2) + 1)>1e-6 & abs(original_data(:,3) + 1)>1e-6 & original_data(:,4)>0& data.drug_redfern(:) > 1);
else
    use_this_data = find(abs(original_data(:,1) + 1)>1e-6 & abs(original_data(:,2) + 1)>1e-6 & abs(original_data(:,3) + 1)>1e-6 & original_data(:,4)>0);
end

if combine_cats_1_and_2
    idxes = find(data.drug_redfern(:)==1 | data.drug_redfern(:)==2);
    output_groups(idxes) = 2;
    data.drug_redfern(idxes)=2;
end

if remove_ranolazine
    use_this_data = use_this_data(use_this_data~=12);
end

for i=1:length(use_this_data)
    fprintf('%i \t %s \t %i\n',use_this_data(i),data.drug_names{use_this_data(i)}, data.drug_redfern(use_this_data(i)))
end
pause
original_use_this_data = use_this_data;
% Stratify the data - done alphabetically
if remove_cat_1==false
    strats = {[11 16 33 47 52 90 82], ... % group 1
        [3 20 24 35 41 37 54 71 94], ... % group 2
        [7 19 25 42 49 86 87]}; % group 3
    if remove_ranolazine
        strats{4} = [14 34 29 57 51 89 74]; % group 4    
    else
        strats{4} = [14 34 29 57 51 89 74 12]; % group 4    
    end
else
    strats = {[1 16  33 47 52 90 82], ... % group 1
        [20 24 35 41 37 54 71 94], ... % group 2
        [19 25 42 49 86 87]}; % group 3
    if remove_ranolazine
        strats{4} = [34 29 57 51 89 74]; % group 4    
    else
        strats{4} = [34 29 57 51 89 74 12]; % group 4    
    end
end

for stratification=1:5
    
    if stratification ==1
        use_this_data = [strats{2} strats{3} strats{4}];
    elseif stratification ==2
        use_this_data = [strats{1} strats{3} strats{4}];
    elseif stratification ==3
        use_this_data = [strats{1} strats{2} strats{4}];
    elseif stratification ==4
        use_this_data = [strats{1} strats{2} strats{3}];
    else
        use_this_data = original_use_this_data;
    end
    
    % Select only this data
    input_data = original_data(use_this_data,:);
    
    % 'Output' - the actual categories
    redfern_cats = data.drug_redfern(use_this_data);
    
    %% Just display the drugs we are going to use with their category...
    fprintf('Running classification for %i drugs:\n',length(use_this_data))
    for i=1:length(use_this_data)
        fprintf('Drug = %s,\t Redfern Category = %i\n',data.drug_names{use_this_data(i)},redfern_cats(i))
    end
    fprintf('\n')
    
    best_mean = 5;
    best_std = 5;
    best_x = 0;
    best_y = 0;
    for x=1:length(original_names)
        x
        for y=x+1:length(original_names)
            %% CHOOSE WHICH INPUTS TO USE FOR DISCRIMINATION
            if dimension_of_analysis == 1
                %disp(['Running ' original_names{x}])
                input_choice_data = [input_data(:,x)];
                labels = {original_names{x}};
            elseif dimension_of_analysis==2
                %disp(['Running with ' original_names{x} ' and ' original_names{y}])
                input_choice_data = [input_data(:,x) input_data(:,y)];
                labels = {original_names{x}, original_names{y}};
            end
            title_text = 'Classification';
            
            % This has to be here in case the operations above reorder the
            % variables.
            use_these_vars = 1:dimension_of_analysis;
            input_choice_data = input_choice_data(:,use_these_vars);
            
            %% Run leave-one-out classification analysis
            error_score = 0;
            %subplot(2,2,1)
            for run=1:length(use_this_data)
                [v_score] = feval(whose_code,run, input_choice_data, redfern_cats,fit_type,plotting_on);
                error_score(run) = v_score;
                %             figure(1)
                %             xlabel(labels{1})
                %             if dimension_of_analysis>1
                %                 ylabel(labels{2})
                %             end
            end
            title(title_text)
            
            % This score is the mean difference between the predicted value and the
            % actual value. It will work for any discrimination technique.
            %fprintf('\nResult: Mean Error: %g , Std Dev: %g \n',mean(error_score),std(error_score))
            if dimension_of_analysis == 1
                mean_error_in_category(x) = mean(error_score);
                std_dev_error_in_category(x) = std(error_score);
            else
                mean_error_in_category(x,y) = mean(error_score);
                std_dev_error_in_category(x,y) = std(error_score);
            end
            
            new_best = false;
            if mean(error_score) < best_mean
                new_best=true;
                best_mean = mean(error_score);
                best_x = x;
                fprintf('Improvement with %s',original_names{x})
                if dimension_of_analysis>1
                    best_y = y;
                    fprintf(' and %s',original_names{y})
                end
                fprintf('\n')
                
                subplot(2,2,2)
                xlabel(labels{1})
                if dimension_of_analysis>1
                    ylabel(labels{2})
                end
                title('Best Predictors')
            end
            
            if new_best
                do_these = [3 4]; % Replot the current graph on the right.
            else
                do_these = 3;
            end
            
            for i=do_these
                subplot(2,2,i)
                hist(error_score)
                xlabel('Error in Category')
                ylabel('Number of drugs')
                ylim([0 23])
                xlim([0 max(redfern_cats)-min(redfern_cats)])
                set(gca,'XTick',0:4)
                title(['Mean = ' num2str(mean(error_score))...
                    ', Std Deviation = ' num2str(std(error_score))])
                
            end
        end
    end
    
    % Save the results of the classifications into structures for messing with
    % in plot_how_good_predictions_are.m
    if stratification <=4
        extra_name = ['_strat_' num2str(stratification) '_'];
    else
        extra_name = '_';
    end
    if remove_cat_1
        extra_name = [extra_name 'no_cat_1'];
    end
    if combine_cats_1_and_2
        extra_name = [extra_name 'combine_12'];
    end
    save(['mean_error_in_category_' num2str(dimension_of_analysis) 'D' extra_name '.mat'],'mean_error_in_category','-mat')
    save(['std_dev_error_in_category_' num2str(dimension_of_analysis) 'D' extra_name '.mat'],'std_dev_error_in_category','-mat')
    
end
