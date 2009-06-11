% Load ClonalPopulations.dat
close all
clear

x = load('ClonalPopulations2.dat');

% Start at time zero
x(:,1) = x(:,1)-300;

createClonalPopulationfigure(x(:,1),x(:,2))