function [ training_points validation_data actual_group ] = convert_training_data_to_cell_groups( data, groups, varargin )
%CONVERT_TRAINING_DATA_TO_CELL_GROUPS Summary of this function goes here
%   Detailed explanation goes here

if nargin==3
    leave_out_point = varargin{1};
else
    leave_out_point = -1;
    assert(nargout==1)
end

% Count how many different groups there are
a = intersect(groups,1:max(groups));
num_groups = length(a);

% Convert the input data into the required format
training_points = cell(num_groups);
for i=1:num_groups
    training_indices=[];
    for j=1:length(groups)
        if (groups(j)==a(i))
            if (j~=leave_out_point)
                training_indices=[training_indices; j];
            else
                validation_indices = j;
                actual_group = groups(j);
            end
        end
        
    end
    training_points{i} = data(training_indices,:);
end

if nargout>=2
    validation_data = data(validation_indices,:);
end