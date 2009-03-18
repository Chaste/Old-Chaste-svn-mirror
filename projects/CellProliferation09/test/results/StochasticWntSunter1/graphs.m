% Calls the plotting function specifying a title and this file path.
close all
clear

title_string = 'Wnt ODE Model Cells in Sunter i) Geometry';
this_path = pwd;

addpath('../');
plot_40_9_data(title_string, this_path);