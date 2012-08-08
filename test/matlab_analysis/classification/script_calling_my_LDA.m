function script_calling_my_LDA( input_args )
%SCRIPT_CALLING_MY_LDA Summary of this function goes here
%   Detailed explanation goes here
close all

% Get some points from different gaussian distributions
% For two of the four variables each class will use the same mean


% Here we choose what approach to take:
% We can recast into linear combinations of the original variables
% Otherwise we just re-order the original variables...
transform_variables = false; 

% We choose whether to generate new points from random normal distributions
new_points = false;

%% Generate some training data...
if new_points    
    num_points_each_class = [8 5 10 12];
    K = length(num_points_each_class);
    base_means = [2 4 2 0];
    class_var = [0 0 0 -1;
        0 5 0 2;
        0 1 0 -1;
        0 9 0 2];
    
    num_vars = length(base_means);
    training_points = cell(K);
    for k = 1:K
        r = zeros(num_points_each_class(k),num_vars);
        for discrim_idx = 1:num_vars
            mean_val = base_means(discrim_idx) + class_var(k,discrim_idx);
            std_val = 1.0;
            r(:,discrim_idx) = mean_val + std_val.*randn(num_points_each_class(k),1)';
        end
        training_points{k} = r;
    end
    save('training_points.mat','training_points','-mat');
else
    load('training_points.mat')
end



%% Cast into the most discriminatory coordinate system...
K = length(training_points); % Number of classes

[transformed_training_points transformation_mat ordering] = ...
    define_reduced_rank_ordering(training_points)

% Choose to use new canonical variables 
if transform_variables
    training_points = transformed_training_points;
else % or use the ordering of the discriminant variables only
    for k=1:K
        this_set = training_points{k};
        this_set = this_set(:,ordering);
        training_points{k} = this_set;
    end
end


%% Choose which of the 4 variables to discriminate on.
% FIRST TWO are always the best after the re-ordering/transform has happened...
use_these_variables = [1 2];

% Select only the variables we are interested in discriminating over.
groups = [];
for i=1:K
    temp = training_points{i};
    training_points{i} = temp(:,use_these_variables);
    groups = [groups; ones(size(temp,1),1)*i];
end



%% Plot the training data
redfern_gscatter_wrapper([],training_points, groups,6)
hold on

x_lims = get(gca,'xlim');
y_lims = get(gca,'ylim');

%% Define the discriminant functions using the training points.
discriminant_functions = define_discriminant_functions(training_points,'linear');

%% Loop over two dimensions of the whole space, make some test points
K = length(training_points);
p = size(training_points{1},2);
resolution=100;
k=0;
test_point = [];
for i=1:resolution
    for j=1:resolution
        
        x(i,j) = x_lims(1)+(x_lims(2)-x_lims(1))*i/resolution;
        y(i,j) = y_lims(1)+(y_lims(2)-y_lims(1))*j/resolution;
        
        point = [x(i,j) y(i,j)];
        test_point = [test_point; point];
    end
end

%% Classify the points
assert(size(test_point,2)==p)
a = classify_these_points(discriminant_functions, test_point );

%% plot them on the graph
k=0;
for i=1:resolution
    for j=1:resolution
        k=k+1;
        % Plot this point on the scatter graph.
        if a(k)==1
            marker = 'k.';
        elseif a(k)==2
            marker = 'r.';
        elseif a(k)==3
            marker = 'm.'; 
        else
            marker = 'b.';
        end
        plot(x(i,j),y(i,j),marker,'MarkerSize',3)
    end
end

%% Plot the location of the class centroids.
M = evaluateCentroids(training_points); 
for i=1:K
    if i==1
        marker = 'ko';
    elseif i==2
        marker = 'ro';
    elseif i==3
        marker = 'mo';
    else
        marker = 'bo';
    end
    plot(M(i,1), M(i,2), marker,'MarkerSize',12);
end

function M = evaluateCentroids(training_points)
% This method is just for plotting not calculation of anything useful
p = size(training_points{1},2); % # of discriminant vars
for k=1:length(training_points)
    temp = training_points{k};
    for j=1:p
        M(k,j) = mean(temp(:,j));
    end
end
