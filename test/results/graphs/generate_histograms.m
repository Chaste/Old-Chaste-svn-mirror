function [ output_args ] = generate_histograms( input_args )
%GENERATE_HISTOGRAMS Generate 3D histograms of crypt invasion data
%   For each mutation state, considering only mutations occurring in the
%   bottom 1/20th of the crypt and not affecting cell viscosity, plot a 3D
%   histogram of the frequency of 'successful' (at taking over the
%   crypt) mutations at given initial heights taking given times. Not
%   necessarily science, but it does show how to use hist3.
data = loadData();
% Example dataset
figure
for i=1:4
    subplot(2,2,i)
    dataset = data{i,1,3};

    takeover_index = find(dataset(:,4)==1);
    takeover_data = dataset(takeover_index,2:3)

    takeover_data(:,1) = takeover_data(:,1)*100/22; % convert to percentage
    takeover_data(:,2) = takeover_data(:,2)/24; % convert to days
    hist3(takeover_data), axis square
    xlabel('Initial mutation height (% of crypt height)')
    ylabel('Experiment duration (days)')
    zlabel('Frequency')
    title(['Mutation state = ' num2str(i)])
end