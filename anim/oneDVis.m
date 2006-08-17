% Generate 1D dataset
close all
clear

load data.dat
[NumberOfTimeSteps NumPoints]=size(data);
NumPoints = NumPoints-1;

for i=1:NumberOfTimeSteps
    plot(data(i,2:end),zeros(NumPoints))
    hold on
    plot(data(i,2:end),zeros(NumPoints),'x')
    axis([0 1.3 -0.1 0.1])
    hold off
    M(i) = getframe;
end

movie(M,10)
