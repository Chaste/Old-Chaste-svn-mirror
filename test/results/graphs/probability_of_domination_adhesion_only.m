function [ this_bucket_height success_prob errors ] = probability_of_domination_adhesion_only( input_args )
%PROBABILITY_OF_DOMINATION Summary of this function goes here
%   Detailed explanation goes here
[data, settings] = loadData();

for i=1:length(settings.viscosity_names)
    data_reduced{i} = data{4,1,i};
end

% Put the height into sub-boxes...
num_viscos = length(settings.viscosity_names);

% This is option for how many sub-divisions to use
num_buckets = 1; % Per proliferation height (anything else doesn't make sense)

% These are hard-coded into the simulations
crypt_length = 14.500;

this_bucket_height=[];
for i=1:num_viscos
    this_bucket_height(i) = str2num(settings.viscosity_names{i});
end
this_bucket_height = this_bucket_height';

num_bucket = zeros(num_buckets*num_viscos,1);
success_bucket = zeros(num_buckets*num_viscos,1);
failure_bucket = zeros(num_buckets*num_viscos,1);

for viscosity_idx = 1:length(settings.viscosity_names)
    try
        dataset = data_reduced{viscosity_idx};
    catch
        break;
    end
    if ~isempty(dataset)

        for j=1:num_buckets
            i = j + (viscosity_idx-1)*num_buckets; % global index
            this_bucket_dataset = dataset;
            num_bucket(i) = num_bucket(i) + length(dataset);
            success_bucket(i) = success_bucket(i) + length(find(this_bucket_dataset(:,4)==1));
            failure_bucket(i) = failure_bucket(i) + length(find(this_bucket_dataset(:,4)==0));
            % In the case of 'no result' we want to remove these
            % experiments...
            idx = find(this_bucket_dataset(:,4)<-0.5);
            if (~isempty(idx))
                num_bucket(i)=success_bucket(i)+failure_bucket(i);
            else            
                assert(num_bucket(i)==success_bucket(i)+failure_bucket(i));
            end
        end
    end
end

assert(length(num_bucket)==length(settings.viscosity_names));

% Plot of frequencies for different heights
figure(1)
bar(this_bucket_height,[success_bucket failure_bucket],'Stacked')
xlabel('Mutant adhesion parameter \alpha','FontSize',16)
ylabel('Frequency','FontSize',16)
hold on
legend('Success','Failure')
set(gca,'FontSize',14)

% Plot of probabilities for different heights
figure(2)
success_prob = (success_bucket./num_bucket);
failure_prob = (failure_bucket./num_bucket);

fprintf('Viscosity\tNumber of Experiments\n');
for i=1:length(num_bucket)
    fprintf('%s \t %g \n',settings.viscosity_names{i},num_bucket(i));
end

confidence_interval = 0.975; % 95% CI altogether for both sides
% for standard normal (central limit thm)
errors = norminv(confidence_interval,0,1) .* sqrt((success_prob.*(1-success_prob))./num_bucket);

success_prob = 100.0*success_prob;
hold on
errorbar(this_bucket_height,success_prob,100*errors,'kx')
xlabel('Mutant adhesion parameter \alpha','FontSize',16)
ylabel('Probability (%)','FontSize',16)
set(gca,'FontSize',14)




