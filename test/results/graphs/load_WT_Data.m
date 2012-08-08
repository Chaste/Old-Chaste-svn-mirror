function [data]= load_WT_Data( input_args )
%LOAD_WT_DATA Summary of this function goes here

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

base_filename = '../processed/processed_results_WT';

full_name = [base_filename '.dat'];

dataset = load(full_name);

data = reshape(dataset, [3 length(dataset)/3])';
disp([full_name ' loaded']);
