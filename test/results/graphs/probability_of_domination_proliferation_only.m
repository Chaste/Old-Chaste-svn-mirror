function [ Z ] = probability_of_domination_proliferation_only( input_args )
%PROBABILITY_OF_DOMINATION Summary of this function goes here

% Get the 'control' simulation from old results
data = loadData();
control_data = data{4,1,7};
control_prob = sum(control_data(:,4))./length(control_data);

%   Detailed explanation goes here
[data, settings] = load_proliferation_only_data();

% Put the height into sub-boxes...
num_heights = length(settings.proliferation_heights);

% These are hard-coded into the simulations
crypt_length = 14.500;

a = str2num(settings.proliferation_heights{1});
b = str2num(settings.proliferation_heights{end});
c = num_heights;
this_bucket_height = linspace(a,b,c)';


num_bucket = zeros(num_heights,1);
success_bucket = zeros(num_heights,1);
failure_bucket = zeros(num_heights,1);

difference_in_force_dominate = zeros(num_heights,1);
difference_in_force_washed = zeros(num_heights,1);
difference_in_force = zeros(num_heights,1);

average_mutant_force_dominate = zeros(num_heights,1);
average_mutant_force_washed = zeros(num_heights,1);
average_wt_force_dominate = zeros(num_heights,1);
average_wt_force_washed = zeros(num_heights,1);
average_mutant_force = zeros(num_heights,1);
average_wt_force = zeros(num_heights,1);

for proliferation_height = 1:length(settings.proliferation_heights)
    try
        dataset = data{proliferation_height};
    catch
        break;
    end
    if ~isempty(dataset)
        i = 1 + (proliferation_height-1); % global index

        num_bucket(i) = num_bucket(i) + length(dataset);
        
        success_idxs = find(dataset(:,8)==1);
        fail_idxs = find(dataset(:,8)==0);   
        
        % % Or look at the simulations where there was a decent fight.
        %success_idxs = find(dataset(:,8)==1 & dataset(:,7)>500);
        %fail_idxs = find(dataset(:,8)==0 & dataset(:,7)>500);   
        
        % Average forces were already printed "per unit time" on bottom of
        % crypt, so here multiply by that time, sum and then divide by
        % total time.
        
        difference_in_force_dominate(i) = sum(-(dataset(success_idxs,6) - dataset(success_idxs,5)).*dataset(success_idxs,7))./sum(dataset(success_idxs,7));
        difference_in_force_washed(i) = sum(-(dataset(fail_idxs,6) - dataset(fail_idxs,5)).*dataset(fail_idxs,7))./sum(dataset(fail_idxs,7));
        difference_in_force(i) = sum(-(dataset(:,6) - dataset(:,5)).*dataset(:,7))./sum(dataset(:,7));

        average_mutant_force_dominate(i) = sum(dataset(success_idxs,6).*dataset(success_idxs,7))./sum(dataset(success_idxs,7));
        average_mutant_force_washed(i) = sum(dataset(fail_idxs,6).*dataset(fail_idxs,7))./sum(dataset(fail_idxs,7));
        average_wt_force_dominate(i) = sum(dataset(success_idxs,5).*dataset(success_idxs,7))./sum(dataset(success_idxs,7));
        average_wt_force_washed(i) = sum(dataset(fail_idxs,5).*dataset(fail_idxs,7))./sum(dataset(fail_idxs,7));
        
        average_mutant_force(i) = sum(dataset(:,6).*dataset(:,7))./sum(dataset(:,7));
        average_wt_force(i) = sum(dataset(:,5).*dataset(:,7))./sum(dataset(:,7));

        success_bucket(i) = length(success_idxs);
        failure_bucket(i) = length(fail_idxs);
        %assert(num_bucket(i)==success_bucket(i)+failure_bucket(i));        
        
    end
    clear dataset
end

% Plot of frequencies for different heights
figure(1)
bar(100*this_bucket_height,[success_bucket failure_bucket],'Stacked')
xlabel('Proliferation ceiling y (% of crypt)','FontSize',16)
xlim(100*[-0.1+this_bucket_height(1) this_bucket_height(end)+0.1])
ylabel('Frequency','FontSize',16)
hold on
legend('Success','Failure')
set(gca,'FontSize',14)

% Plot of probabilities for different heights
figure(2)
num_bucket
success_prob = (success_bucket./num_bucket);
% Hack in the control case
success_prob = [success_prob(1); control_prob; success_prob(2:end)];

confidence_interval = 0.975; % 95% CI altogether for both sides
% for standard normal (central limit thm)
errors = norminv(confidence_interval,0,1) .* sqrt((success_prob.*(1-success_prob))./[num_bucket(1); length(control_data); num_bucket(2:end)]);

hold on
% Add in control
errorbar(100*[this_bucket_height(1); 0.35; this_bucket_height(2:end)],100.*success_prob,100.*errors,'kx--')
xlabel('Proliferation ceiling{\it y}_{thr} (% of height of crypt)','FontSize',16)
xlim(100*[-0.1+this_bucket_height(1) this_bucket_height(end)+0.1])
ylim(100*[0 0.2])
ylabel('Probability of domination (%)','FontSize',16)
set(gca,'FontSize',14)

figure(3)
subplot(1,2,1)
plot(100*this_bucket_height,average_mutant_force_washed,'ko--',...
     100*this_bucket_height,average_wt_force_washed,'ko-')
xlabel('Proliferation ceiling{\it y}_{thr}  (% of height of crypt)','FontSize',16)
ylabel('Average force on cells at base','FontSize',16)
ylim([-4.2 -2.9])
xlim([30 100])
title('Mutant population lost','FontSize',16)
legend('Mutant','WT','Location','SouthWest')
set(gca,'FontSize',14)

subplot(1,2,2)
plot(100*this_bucket_height,average_mutant_force_dominate,'ko--',...
     100*this_bucket_height,average_wt_force_dominate,'ko-')
xlabel('Proliferation ceiling{\it y}_{thr}  (% of height of crypt)','FontSize',16)
ylabel('Average force on cells at base','FontSize',16)
title('Mutant population dominates','FontSize',16)
ylim([-4.2 -2.9])
xlim([30 100])
%legend('Mutant','WT')
set(gca,'FontSize',14)

figure(4)
plot(100*this_bucket_height,average_mutant_force_washed - average_wt_force_washed,'r-')
xlabel('Proliferation ceiling{\it y}_{thr}  (% of height of crypt)','FontSize',16)
ylabel('Average force on cells at bottom when mutant is washed out','FontSize',16)
set(gca,'FontSize',14)

figure(6)
plot(100*this_bucket_height,average_mutant_force,'r-',...
     100*this_bucket_height,average_wt_force,'b-')
xlabel('Proliferation ceiling{\it y}_{thr}  (% of height of crypt)','FontSize',16)
ylabel('Average force on cells at bottom (all outcomes)','FontSize',16)
legend('Mutant','WT')
set(gca,'FontSize',14)

figure(7)
plot(100*this_bucket_height,difference_in_force_dominate,'r-',...
     100*this_bucket_height,difference_in_force_washed,'b-')
 
figure(8)
plot(100*this_bucket_height,average_mutant_force_dominate,'r-',...
     100*this_bucket_height,average_mutant_force_washed,'b-')
xlabel('Proliferation ceiling{\it y}_{thr}  (% of height of crypt)','FontSize',16)
legend('Mutant dominate','Mutant washed')
set(gca,'FontSize',14)

figure(9)
plot(100*this_bucket_height,average_wt_force_dominate,'r-',...
     100*this_bucket_height,average_wt_force_washed,'b-')
xlabel('Proliferation ceiling{\it y}_{thr}  (% of height of crypt)','FontSize',16)
legend('WT dominate','WT washed')
set(gca,'FontSize',14)

