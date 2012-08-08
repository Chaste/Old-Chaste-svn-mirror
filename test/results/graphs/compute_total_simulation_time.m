function [ output_args ] = compute_total_simulation_time( input_args )
%COMPUTE_TOTAL_SIMULATION_TIME Compute the total simulation time taken
%   Add up the simulation time taken for all crypt invasion simulations

% Load all the data
[data settings] = loadData();

total_time = 0.0;

for state=1:length(settings.mutations)
    for height = 1:length(settings.heights)
        for visc = 1:length(settings.viscosity_names)
            dataset = data{state,height,visc};
            if ~isempty(dataset)
                total_time = total_time + sum(dataset(:,3));
            end
        end
    end
end

% Add Proliferation only data
[prolif_data settings] = load_proliferation_only_data();
for i=1:length(settings.proliferation_heights)
    temp = prolif_data{i};
    total_time = total_time + sum(temp(:,3));
end

% Add wild type data.
WT_data = load_WT_Data();
total_time = total_time + sum(WT_data(:,2));

total_time_in_hours = total_time

total_time_in_years = total_time/(24*365)