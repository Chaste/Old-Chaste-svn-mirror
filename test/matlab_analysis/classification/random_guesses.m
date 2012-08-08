
%close all
clear

combine_cats_1_and_2 = true;
remove_ranolazine = true;
N = 100000;

cd ..
data = load_drug_data();
cd classification

%% Load up all possible inputs
load('compiled_results.mat');
original_data = compiled_results.measures;
original_names = compiled_results.measure_names';
assert(length(original_data)==length(original_names));

%% All three IC50 values and dose data
use_this_data = find(abs(original_data(:,1) + 1)>1e-6 & abs(original_data(:,2) + 1)>1e-6 & abs(original_data(:,3) + 1)>1e-6 & original_data(:,4)>0);

if remove_ranolazine
    use_this_data = use_this_data(use_this_data~=12);
end

if combine_cats_1_and_2
    idxes = find(data.drug_redfern(:)==1 | data.drug_redfern(:)==2);
    output_groups(idxes) = 2;
    data.drug_redfern(idxes)=2;
end

% 'Output' - the actual categories
redfern_cats = data.drug_redfern(use_this_data);
random_guess = randi(4,length(redfern_cats),N)+1;

mean_error = zeros(N,1);
std_error = zeros(N,1);
total_errors = zeros(4,1);

for i=1:N
% Simulate some guesses...
    errors = abs(redfern_cats - random_guess(:,i));
    i
    for j=0:1:3
        total_errors(j+1) = total_errors(j+1) + length(find(errors==j));
    end
    mean_error(i) = mean(errors);
    std_error(i) = std(errors);
end
figure
bar([0 1 2 3],(total_errors)./N,0.5,'k')
xlabel('Error in categorization','FontSize',18)
ylabel('Number of drugs','FontSize',18)
ylim([0 23])
xlim([-0.5 max(redfern_cats)-min(redfern_cats)+0.5])
set(gca,'XTick',0:4,'FontSize',16)
title(['mean = ' num2str(mean(mean_error),'% 4.3f')...
    ', std deviation = ' num2str(mean(std_error),'% 4.3f')],'FontSize',18)


% Make sure we only put points on the scatter plot that haven't already
% been plotted (to stop us having eps files with a million points on them.)
[things_to_plotx things_to_ploty] = remove_repeated_points(mean_error, std_error);
length(things_to_plotx)
figure
scatter(things_to_plotx, things_to_ploty,'k.')
hold on
sim_random_mean = mean(mean_error)
sim_random_std = mean(std_error)

plot([sim_random_mean sim_random_mean],[0 3],'k-')
plot([0 5],[sim_random_std sim_random_std],'k-')
meanDecentMarker = 0.323;
stdDecentMarker = 0.541;
plot([meanDecentMarker meanDecentMarker],[0 3],'k:')
plot([0 5],[stdDecentMarker stdDecentMarker],'k:')
set(gca,'FontSize',14);
xlim([0 2.5])
ylim([0 1.5])
xlabel('Mean category prediction error','FontSize',16)
ylabel('Standard deviation of category prediction error','FontSize',16)

hold on




