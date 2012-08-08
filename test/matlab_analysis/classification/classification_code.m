function [validation_scores] = classification_code( leave_out_point, data, grouping, fit_type, plotting_on,varargin )
%PLAYWITHCLASSIFICATION Summary of this function goes here
%   Detailed explanation goes here
% Options

% If we are doing redfern categories use a special gscatter wrapper to 
% ensure that they are given consistent colours and symbols.
redfern_categories = true;
if nargin==6
    redfern_categories = varargin{1};
end

n = length(data);
group_names = sort(grpstats(data,grouping,'gname')); % For use in figure legends

% Sort according to group (just for comparison with my classification code)
[grouping idx] = sort(grouping);
data=data(idx,:);

% Work out what is training data and validation point.
training_indices=[];
for i=1:length(data)
    if i~=leave_out_point
        training_indices = [training_indices; i];
    end
end

training_x = data(training_indices,1);
training_y = data(training_indices,2);
training_group = grouping(training_indices);
validation_x = data(leave_out_point,1);
validation_y = data(leave_out_point,2);
validation_group = grouping(leave_out_point);

if plotting_on
    % Plot the training data in a scatter plot (big points).
    hold off
    h1 = redfern_gscatter_wrapper(training_x,training_y,training_group,8^2,[]);
    set(h1,'LineWidth',2)
    legend(group_names,'Location','NorthEastOutside')
    hold on
end

% Estimate the classification for validation data.
[C] = classify([validation_x validation_y],[training_x training_y],...
    training_group,fit_type,'empirical'); % Prior depends on freq of training classes
if plotting_on
    % Plot the validation data in the same plot
    h1 = redfern_gscatter_wrapper(validation_x,validation_y,C,12^2,'o');
    set(h1,'LineWidth',2)
end

% Work out the log likelihoods for training and test sets
% Estimates for the actual training dataset
[training_log_likelihood train_num_right train_num_wrong] = ...
    calculate_log_likelihood(training_x,training_y,training_group,...
    training_x,training_y,training_group,fit_type, redfern_categories);
training_log_likelihood = training_log_likelihood/length(training_group);
% Estimates for the validation data
[validation_log_likelihood validation_num_right validation_num_wrong validation_scores] = ...
    calculate_log_likelihood(validation_x,validation_y,validation_group,...
    training_x,training_y,training_group,fit_type, redfern_categories);
validation_log_likelihood = validation_log_likelihood/length(validation_group);


fprintf('Point ');
for i=1:length(data(leave_out_point,:))
    fprintf('%d, ',data(leave_out_point,(i)));
end
fprintf(': Actually %i, classified as %i, error %i \n',validation_group,C, validation_scores);


if plotting_on
    % Estimate the classification for all points (for graphing only)
    min_x = min([training_x; validation_x]);
    max_x = max([training_x; validation_x]);
    min_y = min([training_y; validation_y]);
    max_y = max([training_y; validation_y]);
    
    [Xs,Ys] = meshgrid(linspace(min_x,max_x),...
        linspace(min_y,max_y));
    Xs = Xs(:); Ys = Ys(:);
    [Cs] = classify([Xs Ys],[training_x training_y],training_group,fit_type,'empirical');
    
    redfern_gscatter_wrapper(Xs,Ys,Cs,1,'.');
    xlim([min_x max_x])
    ylim([min_y max_y])
    
    % Plot the decision boundaries...
    % plot_decision_boundaries(coeff,fit_type);
end

%% Calculate the log likelihood of getting the correct classification of
% each data point
function [log_likelihood num_right num_wrong sum_error_scores] = ...
    calculate_log_likelihood(x,y,actual_groups,training_x,training_y,...
                             training_groups,fit_type, redfern_categories)

log_likelihood = 0;

% Do the classification again - same training data.
[C,err,P,logp,coeff] = classify([x y],[training_x training_y],training_groups,fit_type,'empirical');

for i=1:length(actual_groups)
    % Probability of being in the right category =
    % P(i,index_actual_groups(i));
    log_likelihood = log_likelihood + log(P(i,actual_groups(i)-min(actual_groups)+1));
end

num_wrong = 0;
num_right = 0;
sum_error_scores = 0;

for i=1:length(actual_groups)
    if C(i)==actual_groups(i)
        num_right = num_right + 1;
    else
        num_wrong = num_wrong + 1;
    end
    % Here we provide a score of error
    % Get 0 if the categories are the same as predicted
    % 1 if they are 1 away , 2 if 2 away etc.
    sum_error_scores = sum_error_scores + abs(C(i)-actual_groups(i));
end



%% Plot the lines which form the decision boundaries.
function plot_decision_boundaries(coeff, type)
axes_limits = [get(gca,'XLim') get(gca,'YLim')];

% Plot the decision boundary lines.
for i=1:size(coeff,1)-1
    for j=i+1:size(coeff,2)
        K = coeff(i,j).const;
        L = coeff(i,j).linear;

        if strcmp(type, 'quadratic')==1
            Q = coeff(i,j).quadratic;
            f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',...
                K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));

        elseif strcmp(type,'linear')
            f = sprintf('0 = %g+%g*x+%g*y',K,L);
        else
            error('Coefficient type unrecognised')
        end

        figure_handle = ezplot(f,axes_limits);
        set(figure_handle,'Color','m','LineWidth',2)
    end
end