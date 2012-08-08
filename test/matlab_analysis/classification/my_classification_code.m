function [ error_score] = my_classification_code( leave_out_point, data, groups, fit_type,plotting_on )
% MY_CLASSIFICATION_CODE Summary of this function goes here
%   LEAVE_OUT_POINT - allows us to use one point for validation.
%   DATA            - The dataset, each row is one training point 
%                     in class given by row of 'groups';
%   GROUPS          - any integers.

% Sort according to group (just for comparison with classify() code)
[groups idx] = sort(groups);
data=data(idx,:);

[training_points validation_data actual_group] = ...
    convert_training_data_to_cell_groups(data, groups, leave_out_point);

% Count how many different groups there are
% Groups are given by integers
a = intersect(groups,1:max(groups));
num_groups = length(training_points);

for i=1:num_groups
    group_to_redfern_mapping(i) = a(i);
    group_names{i}=num2str(a(i));
end

%%
%%% WE NOW DO ALL THE FANCY RANK REDUCTION BEFORE CALLING THE LEAVE ONE OUT
%%% LOOP
% [transformed_training_points transformation_mat ordering] = ...
%     define_reduced_rank_ordering(training_points)
% 
% If we want to use the original variables in order of discriminant
% % ability:
% disp(['Taking variables ' num2str(ordering(1)) ' & ' num2str(ordering(2))...
%     ' to use as discriminant variables'])
% data(validation_indices,:) = data(validation_indices,ordering);
% for k=1:length(training_points)
%     this_set = training_points{k};
%     this_set = this_set(:,ordering); 
%     training_points{k} = this_set;
% end
% If we want to use the transformed variables
% data = transform_coordinates(data, transformation_mat);
% training_points = transformed_training_points;
% %% Make training and validation data two dimensional
% use_these_vars = [1 2]; % columns 1 and 2 should in theory now be the most useful
% validation_data = data(validation_indices,use_these_vars);
% for k=1:length(training_points)
%     this_set = training_points{k};
%     this_set = this_set(:,use_these_vars); % Get only primary 2 discriminant vars
%     training_points{k} = this_set;
% end

%% Plot the training data
if plotting_on
    redfern_gscatter_wrapper([],training_points,[],64,[],group_to_redfern_mapping);
    legend(group_names,'Location','NorthEastOutside')
    hold on
end

% Define the discriminant functions using the training points.
discriminant_functions = define_discriminant_functions(training_points,fit_type);

% Evaluate the classification of the validation point
[validation_classified_as] = classify_these_points(discriminant_functions,...
    validation_data,group_to_redfern_mapping);

%Here we provide a score of error
%Get 0 if the categories are the same as predicted
%1 if they are 1 away , 2 if 2 away etc.
error_score = abs(validation_classified_as - actual_group);

% Useful to check that the error classification scoring is working
% properly.
% fprintf('Point ');
% for i=1:length(validation_data)
%     fprintf('%d, ',validation_data(i));
% end
fprintf(': Actually %i, classified as %i, error %i \n',actual_group,validation_classified_as, error_score);

if plotting_on
    % Plot the validation data in the same plot (big points)
    redfern_gscatter_wrapper(validation_data(1),validation_data(2)...
        ,validation_classified_as,144,'o');
    
    %% Estimate the classification for all points (for graphing only)
    [Xs,Ys] = meshgrid(linspace(min(get(gca,'Xlim')), max(get(gca,'Xlim')), 100),...
        linspace(min(get(gca,'Ylim')), max(get(gca,'Ylim')), 100));
    Xs = Xs(:); Ys = Ys(:);
    Cs =classify_these_points(discriminant_functions,[Xs Ys],group_to_redfern_mapping);    
    redfern_gscatter_wrapper(Xs,Ys,Cs,1,'.');
    hold off
    
end







