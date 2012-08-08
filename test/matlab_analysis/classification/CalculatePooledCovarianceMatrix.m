function [ sigma_mat pooled_sigma ] = CalculatePooledCovarianceMatrix( training_points )
%CALCULATEPOOLEDCOVARIANCEMATRIX Summary of this function goes here
%   Detailed explanation goes here

K = length(training_points);
dim_measures = size(training_points{1},2); % This is 'p' in the textbook

% Work out how much data we have in each category (NK) and in total (N)
% and scatter plot the training data points

NK=zeros(1,K);

% Set up an empty matrix for pooled_sigma (for LDA).
pooled_sigma = zeros(dim_measures, dim_measures);
sigma_mat = cell(K);

% Loop over each class, work out mean and covariance matrix
for i=1:K
    NK(i) = size(training_points{i},1);
    sigma_mat{i} = CalculateCovarianceMatrix(training_points{i});
    % Add to pooled sigma for LDA scale according to size properly
    pooled_sigma = pooled_sigma + sigma_mat{i}.*(NK(i)-1);
end

% This equation according to calculations below (4.10)
pooled_sigma = pooled_sigma./(sum(NK) - K);



