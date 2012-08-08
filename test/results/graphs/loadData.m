function [data settings]= loadData( input_args )
%LOADDATA Summary of this function goes here
%   The data is returned as a cell array for
% data{m,h,v}
% m = mutation state 1/2/3/4
% h = height 0/1/2
% v = viscosity 0.5/0.6/0.7/0.75/0.8/0.9/1.0/1.25/etc. see below
% 
%   The data in each cell has as its columns
% 1. Experiment number (unsigned)
% 2. Height of initial mutation (double)
% 3. Duration (double)
% 4. Success? (bool)
% 
close all

base_filename = '../processed/processed_results';

settings.mutations = {'1' '2' '3' '4'};
settings.mutation_names = {'Proliferation 50%' 'Proliferation 100%' 'Proliferation 90%' 'Proliferation 35% (wild type)'};
settings.proliferation_heights = {'0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0'};
settings.heights = {'0' '1' '2'};
settings.viscosity_names = {'0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0' '1.25' '1.5' '1.75' '2.0' '2.25' '2.5' '2.75' '3.0' '3.25' '3.5' '3.75' '4.0' '4.5' '5.0' '5.5' '6.0' '7.0' '8.0' '9.0' '10.0'};
for i=1:length(settings.viscosity_names)
    settings.viscosities(i) = str2double(settings.viscosity_names{i});
end

for m = 1:length(settings.mutations)
    for h = 1:length(settings.heights)
        for v = 1:length(settings.viscosity_names)
            try
                full_name = [base_filename settings.mutations{m} '_' ...
                    settings.heights{h} '_' settings.viscosity_names{v} '.dat'];
                dataset = load(full_name);
                disp([full_name ' loaded'])
                data{m,h,v} = reshape(dataset, [4 length(dataset)/4])';
            catch
                disp([full_name ' doesn''t exist']) 
                data{m,h,v} = [];
            end            
        end
    end
end
