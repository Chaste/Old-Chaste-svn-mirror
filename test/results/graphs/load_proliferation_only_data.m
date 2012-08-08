function [data settings]= load_proliferation_only_data( input_args )
%LOAD_PROLIFERATION_ONLY_DATA Summary of this function goes here
%   The data is returned as a cell array for
% data{h}
% h = proliferation height 0.3/0.4/0.5/0.6/0.7/0.8/0.9/1.0
% 
%   The data in each cell has as its columns
% 1. Experiment number (unsigned)
% 2. Height of initial mutation (double)
% 3. Duration (double)
% 4. Number of cells in crypt at end (unsigned)
% 5. Average force on WT
% 6. Average force on mutant
% 7. Time which was averaged over (time mutants at base)
% 8. Domination? (bool)
% 
close all

base_filename = '../processed/processed_results_proliferation_only';

settings.proliferation_heights = {'0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0'};
settings.heights = {'0'};
for h = 1:length(settings.proliferation_heights)
    full_name = [base_filename settings.proliferation_heights{h} '.dat'];
    try        
        dataset = load(full_name);    
    catch
        disp([full_name ' doesn''t exist']) ;
    end
    try
        data{h} = reshape(dataset, [8 length(dataset)/8])';
    catch
        disp([full_name ' appears to be incomplete.']);
    end
    disp([full_name ' loaded']);
end
