function [U D invD] = eigenvalue_decomposition(sigma)

[U D] = schur(sigma); % Inbuilt matlab eigenvalue decomposition method

% Get inverse of diagonal matrix whilst we're here.
invD = zeros(length(D),length(D));
for d=1:length(D)
    invD(d,d) = 1.0./D(d,d);
end