%
% Code to plot the output of the Chaste test
% TestMeinekeLabellingExperiments.hpp
%
% Requires files of the number of cells in each position that were
% labelled in a column, at time = 40 minutes and time = 9 hours.
% 
% Author: Gary Mirams
% Date: 9/11/07

close all
clear

title_string = 'Meineke Cells in Sunter iii) Geometry';
% This is the number of separate runs of the program
runs = 4;
% This is the number of simulations each run did.
simulations = 51;

% Don't mess below 'ere!
sections_taken = runs*simulations*2;

t40_results = zeros(100,1);
t9_results = zeros(100,1);

for i=1:runs

    sim_str = num2str(i);
    simulations_str = num2str(simulations);

    t40_results = t40_results + load(['labelled_t40_totals', simulations_str, '_run_', sim_str, '.dat']);
    t9_results = t9_results + load(['labelled_t9_totals', simulations_str, '_run_', sim_str, '.dat']);

end

% Scale to percentages
t40_results = 100.*t40_results./sections_taken;
t9_results = 100.*t9_results./sections_taken;

plot(t40_results,'r.')
hold on
plot(t9_results,'bx')
xlabel('Cell Position')
ylabel('% Labelled Cells')
title([title_string , ': ' , num2str(runs*simulations), ' experiments'])
legend('t = 0:40','t = 9:00')
xlim([0 50])
ylim([0 100])

