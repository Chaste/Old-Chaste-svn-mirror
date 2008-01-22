% PlotCellPerimsAndAreas.m
% It reads in the positions of all cells at the beginning and end of a 
% Meineke-stlye labelling experiment and plots the percentages of cells 
% that are labelled in ranges 
%

close all
clear

% Experiment Setup
title_string = 'Simple Wnt Cells in Sunter i) Geometry';
crypt_height = 30;
path = '/tmp/pmxaw/testoutput/Noddy_WNT_Yes_Area_Yes_Length/results_from_time_4/vis_results/';
% End of setup

addpath('../');	% Adds the LoadNonConstantLengthData function.

Voronoi_data = LoadNonConstantLengthData([path 'results.visvoronoi']);

X = [];
Y = [];
Area = [];
Perimeters = [];




for i=1:length(Voronoi_data) % time loop
	
	num_cells = (length(Voronoi_data{i})-1)/5;
	for j = 1:num_cells
		X = [X Voronoi_data{i}(5*j-2)];
		Y = [Y Voronoi_data{i}(5*j-1)];
		Area = [Area Voronoi_data{i}(5*j)];
		Perimeters = [Perimeters Voronoi_data{i}(5*j+1)];
		
	end
	
	
end

figure;
BoxPlotChaste(Area);
figure;
BoxPlotChaste(Perimeters);
