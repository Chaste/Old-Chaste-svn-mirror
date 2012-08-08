function [ Z vecs ordering] = define_reduced_rank_ordering( training_points )
%DEFINE_REDUCED_RANK_ORDERING Get optimal classification directions
%
% Taken from section 4.3.3 of Hastie, Tibshirani & Friedman 
% The Elements of Statistical Learning.
% 
%
% training_points - a cell array of the different classes (length K)
%                   Entries are matrices with:
%                    - width p = Number of discriminant variables.
%                    - length Nk = Number of points in this class.
%
% Outputs
% Ms = Transformed centroid locations
% Z = Transformed training_points into new variables.
% vecs = Transformation matrix (first column is fist canonical variable etc.)
% ordering = as an alternative to using linear combinations to transform
%               this is a ranking of the best original variables to use.
%
% Author: Gary Mirams
% Date: 29/3/2010
%

order = 1:size(training_points{1},2); % Number of discriminant variables

p = size(training_points{1},2); % Number of discriminant variables
K = length(training_points); % Number of classes

order = 1:p; % Default ordering

for k=1:length(training_points)
    temp = training_points{k};
    for j=1:p
        M(k,j) = mean(temp(:,j));
    end
end

[ sigma_mat pooled_sigma ] = CalculatePooledCovarianceMatrix( training_points );

W = pooled_sigma; % I think...

[U D Dinv] = eigenvalue_decomposition(W);

%% Get result via Fisher's method see discussion below (4.16)
B = CalculateCovarianceMatrix(M, true); % Wikipedia says used biased one here.
WinvB = U*Dinv*U'*B;

[vecs] = eigenvalue_decomposition(WinvB);

[q r ordering] = qr(WinvB,0);

% Perform coordinate transformation

Z = cell(K,1);
for k=1:K
    Z{k} = transform_coordinates(training_points{k}, vecs);
end



