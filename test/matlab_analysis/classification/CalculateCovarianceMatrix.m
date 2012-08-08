function [ sigma ] = CalculateCovarianceMatrix( A ,flag )
%CALCULATECOVARIANCEMATRIX Summary of this function goes here
%   Detailed explanation goes here

biased = false;
if nargin==2
    biased=flag;
end

m = size(A,1);
assert(m>1);

% This is just a fast coded matlab method for subtractions.
Ac = bsxfun(@minus,A,sum(A,1)/m);  % Remove mean

% We generally want the unbiased estimate...
if biased
    sigma = (Ac' * Ac) / m;
else
    sigma = (Ac' * Ac) / (m-1);
end