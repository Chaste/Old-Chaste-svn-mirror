function [ Z ] = probability_of_domination( input_args )
%PROBABILITY_OF_DOMINATION Summary of this function goes here
%   Detailed explanation goes here
[data, settings] = loadData();

% Put the height into sub-boxes...
num_mutations = length(settings.mutations);
num_heights = length(settings.heights);
num_viscosities = length(settings.viscosities);

% This is option for how many sub-divisions to use
num_buckets = 10; % Looks reasonably smooth...

% These are hard-coded into the simulations
crypt_length = 14.500;
box_proportion = 1.0/20.0;

big_box_height = box_proportion*crypt_length;
bucket_height = big_box_height/num_buckets;

this_bucket_height_min = linspace(0,big_box_height*num_heights,num_buckets*num_heights)';

for mutation=1:num_mutations
    for viscosity_index = 1:num_viscosities
        num_bucket = zeros(num_buckets*num_heights,1);
        success_bucket = zeros(num_buckets*num_heights,1);
        failure_bucket = zeros(num_buckets*num_heights,1);
        
        for height = 1:num_heights
            start_height = (height-1)*big_box_height;
            try
                dataset = data{mutation,height,viscosity_index};
            catch
                break;
            end
            if ~isempty(dataset)

                for j=1:num_buckets
                    i = j + (height-1)*num_buckets; % global index
                    this_bucket_height_max = this_bucket_height_min(i) + bucket_height;

                    this_bucket_dataset_indices = find(dataset(:,2)>=this_bucket_height_min(i) & ...
                        dataset(:,2)<this_bucket_height_max);
                    this_bucket_dataset = dataset(this_bucket_dataset_indices,:);

                    num_bucket(i) = num_bucket(i) + length(this_bucket_dataset_indices);
                    success_bucket(i) = success_bucket(i) + length(find(this_bucket_dataset(:,4)==1));
                    failure_bucket(i) = failure_bucket(i) + length(find(this_bucket_dataset(:,4)==0));
                    assert(num_bucket(i)==success_bucket(i)+failure_bucket(i));
                end
            end
        end

        % Plot of frequencies for different heights
        figure(1+2*viscosity_index)
        subplot(2,2,mutation)
        bar(this_bucket_height_min,[success_bucket failure_bucket],'Stacked')
        xlabel('Height (of lower bound of bucket)')
        xlim([-0.1 big_box_height*num_heights+0.1])
        title([settings.mutation_names{mutation} ', \alpha = ' settings.viscosity_names{viscosity_index}])
        ylabel('Frequency')
        hold on
        legend('Success','Failure')
        % Plot of probabilities for different heights
        figure(2+2*viscosity_index)
        subplot(2,2,mutation)
        success_prob = (success_bucket./num_bucket);
        failure_prob = (failure_bucket./num_bucket);
        bar(this_bucket_height_min,[success_prob failure_prob],'Stacked')
        xlabel('Height (of lower bound of bucket)')
        xlim([-0.1 big_box_height*num_heights+0.1])
        title([settings.mutation_names{mutation} ', \alpha = ' settings.viscosity_names{viscosity_index}])
        ylabel('Probability')
        hold on
        legend('Success','Failure')

        % Save the success probability
        Z(mutation,viscosity_index,:) = success_prob;
    end
end