function delta_funcs = define_discriminant_functions(training_points, type)
%
% This method defines delta (discriminant) functions for each of the
% K inputted training groups.
%
% It returns K handles to the relevant discriminant functions.
% If you call each discriminant function in turn the one which 
% has the maximum value is the group into which your point will be 
% classified.
%
% From Elements of Statistical Learning, Edition 1, P 87, eqn (4.10).
%
% Author: Gary Mirams
% Date: 20/3/10
%

K = length(training_points);
dim_measures = size(training_points{1},2); % This is 'p' in the textbook

[sigma_mat pooled_sigma] = CalculatePooledCovarianceMatrix(training_points);

% Work out how much data we have in each category (NK) and in total (N)
% and scatter plot the training data points
N=0;
for i=1:K
    NK(i) = size(training_points{i},1);
    N = N + NK(i);
end

for i=1:K
    mean_this_set{i} = mean(training_points{i});
    % Calculate Priors according to number of drugs in each category.
    %pi_this_set{i} = NK(i)./N;    
    % ALTERNATIVELY 
    % let them all be equal proportions of the number of categories
    pi_this_set{i} = 1.0/K;
end

% Work out pooled sigma decomposition just once
[Upool Dpool invDpool] = eigenvalue_decomposition(pooled_sigma);

% Define handles to discriminant functions for each of the classes
for i=1:K
    if strcmp(type,'linear')==1
        % Do some evaluations once here...
        part1(:,i) = invDpool*(Upool'*mean_this_set{i}'); % a vector for each class
        part2(i) = -0.5*(Upool'*mean_this_set{i}')'*invDpool*...
                        (Upool'*mean_this_set{i}'); % a scalar for each class

        % Equation (4.10) using the hints in 
        % section 4.3.2 (Computations for LDA).
        delta_funcs{i} = @(input_data) (...
            (Upool'*input_data')'*part1(:,i) ...
            + part2(i) ...
            + log(pi_this_set{i}));
        
    elseif strcmp(type,'quadratic')==1
        % Equation (4.12) using the hints in 
        % Section 4.3.2 (Computations for LDA).
        [U D invD] = eigenvalue_decomposition(sigma_mat{i});
        
        delta_funcs{i} = @(input_data) (...
            -0.5*(sum(log(diag(D)))) ...
            -0.5*(U'*(input_data' - mean_this_set{i}'))'...
            *invD*(U'*(input_data' - mean_this_set{i}')) ...
            + log(pi_this_set{i}));
    else
        error('Unrecognised discriminant analysis requested.')
    end
end



