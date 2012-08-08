function [ class ] = classify_these_points(discriminant_functions, ...
                                           test_point, ...
                                           group_mapping)
%CLASSIFY_THIS_POINT Use the given discriminant functions to classify a
%point
K = length(discriminant_functions);

for i=1:size(test_point,1)
    point = test_point(i,:);
    
    % Evaluate the discriminant functions for each of the categories.
    discriminant_scores = zeros(1,K);
    for k=1:K
        discriminant_scores(k) = feval(discriminant_functions{k},point);
    end
    
    % Assign this point to the maximum valued discriminant function.
    [val class(i)] = max(discriminant_scores);
end

if nargin==3
    for i=1:length(class)
        class(i) = group_mapping(class(i));
    end
end