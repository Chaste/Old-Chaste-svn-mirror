function [ output_args ] = example_graph( input_args )
%EXAMPLE_GRAPH Summary of this function goes here
%   Detailed explanation goes here
[data, settings] = loadData();
% Example dataset
figure
for i=1:4
    subplot(2,2,i)
    dataset = [data{i,1,3}; data{i,2,3}; data{i,3,3}];

    takeover_index = find(dataset(:,4)==1);
    failure_index = find(dataset(:,4)==0);

    scatter(dataset(failure_index,2), dataset(failure_index,3),'b.')
    hold on
    scatter(dataset(takeover_index,2), dataset(takeover_index,3),'rx')
    xlabel('Mutation height (cell diameters)','FontSize',16)
    xlim([0 2])
    ylabel('Experiment duration (hours)','FontSize',16)
    title(['Mutation state: ' settings.mutation_names{i}],'FontSize',16)
    set(gca,'FontSize',14)
end

%% Export this figure to the paper images directory.
papersize = get(gcf, 'PaperSize');
width = 30.0;         % Initialize a variable for width.
height = 20.0;      % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

print -depsc2 ../../../paper/images/example_duration