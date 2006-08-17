% Generate 1D dataset
close all
clear

NumPoints = 10;
Length = 1.0;
NumberOfTimeSteps = 10;
TimeStep = 0.1;
locations = linspace(0,Length,NumPoints)
for t=1:NumberOfTimeSteps
    time = (t-1)*TimeStep;
    data(t,:) = [time locations];
    locations = locations + 0.02;
    
end
save 'data.dat' data -ascii