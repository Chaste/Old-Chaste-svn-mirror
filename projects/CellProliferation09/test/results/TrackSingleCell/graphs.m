% This code needs single_cell_compile_results.sh to have been run on the relevant set
% of results.
%
% It reads in the positions of the cell and its voronoi information through
% time from the Chaste code TestGeneratePlotsOfASingleCell.hpp
%
% Author: Gary Mirams
% Date: 23/4/2008
%
%
close all
clear

addpath('../../../../../anim/matlab/'); % Adds the LoadNonConstantLengthData function.

all_data = LoadNonConstantLengthData('single_voronoi.dat')';

for i=1:length(all_data)
    if length(all_data{i})>1
        time(i) = all_data{i}(1);
        x(i) = all_data{i}(3);
        y(i) = all_data{i}(4);
        area(i) = all_data{i}(5);
        perimeter(i) = all_data{i}(6);
    else
        time(i) = all_data{i};
    end
end

time = time';
x = x';
y = y';
area = area';
perimeter = perimeter';

figure
plot(time(1:length(area)), area,'k.')
xlabel('Time (hrs)', 'Fontsize', 16)
ylabel('Area of a single cell', 'Fontsize', 16)
set(gca, 'Fontsize', 14)
xlim([550 670])
ylim([0 1])

figure
scatter(y, area, 'k.')
xlabel('Height up crypt (cells)', 'Fontsize', 16)
ylabel('Area of a single cell', 'Fontsize', 16)
set(gca, 'Fontsize', 14)
ylim([0 1])