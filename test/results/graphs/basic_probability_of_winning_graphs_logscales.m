close all
clear all

data = loadData();

% Load just the Wild-Type data (labelled no change to viscosity or
% proliferation) from all heights
data = [data{4,1,7}; data{4,2,7}; data{4,3,7}];

success_indices = find(data(:,4)==1);
failure_indices = find(data(:,4)==0);

num_expts = length(data)

% Sanity check
assert(length(success_indices) + length(failure_indices)==length(data));

num_successes = length(success_indices)

figure
subplot(2,2,2)
scatter(log10(data(success_indices,2)),data(success_indices,3)./24,'k.')
xlim([-7 1])
ylim([0 250])
%xlabel('log_{10} Initial mutation height (cell diameters)','FontSize',16)
%ylabel('Simulation duration (days)','FontSize',16)
set(gca,'FontSize',14)
subplot(2,2,4)
scatter(log10(data(failure_indices,2)),data(failure_indices,3)./24,'k.')
xlim([-7 1])
ylim([0 250])
%xlabel('log_{10} Initial mutation height (cell diameters)','FontSize',16)
%ylabel('Simulation duration (days)','FontSize',16)
set(gca,'FontSize',14)
subplot(2,2,1)
scatter(data(success_indices,2),data(success_indices,3)./24,'k.')
xlim([0 2.2])
ylim([0 250])
%xlabel('Initial mutation height (cell diameters)','FontSize',16)
%ylabel('Simulation duration (days)','FontSize',16)
set(gca,'FontSize',14)
subplot(2,2,3)
scatter(data(failure_indices,2),data(failure_indices,3)./24,'k.')
xlim([0 2.2])
ylim([0 250])
xlabel('Initial mutation height (cell diameters)','FontSize',16)
ylabel('Simulation duration (days)','FontSize',16)
set(gca,'FontSize',14)
% 
% figure
% subplot(2,1,1)
% hist(log10(data(success_indices,2)))
% subplot(2,1,2)
% hist(log10(data(failure_indices,2)))

num_buckets = 15;
heights = linspace(log10(min(data(:,2))),log10(max(data(:,2))),num_buckets)';
success_space_buckets = zeros(num_buckets,1);
fail_space_buckets = zeros(num_buckets,1);
success_durations{num_buckets} = [];
fail_durations{num_buckets} = [];

for i = 1:length(data)
    for b = 1:num_buckets
        if (heights(b) < log10(data(i,2)) && log10(data(i,2)) <= heights(b+1))
            if (data(i,4)==1)
                success_durations{b} = [success_durations{b}; data(i,3)];
                success_space_buckets(b) = success_space_buckets(b) + 1;
            else
                fail_durations{b} = [fail_durations{b}; data(i,3)];
                fail_space_buckets(b) = fail_space_buckets(b) + 1;
            end
        end
    end
end
success_space_buckets
fail_space_buckets
for i=1:num_buckets
    assert(length(success_durations{i})==success_space_buckets(i))
    assert(length(fail_durations{i})==fail_space_buckets(i))
end
probability_space_buckets = success_space_buckets ./ (success_space_buckets + fail_space_buckets);

mean_success_durations  = zeros(1,num_buckets)';
mean_fail_durations = zeros(1,num_buckets)';
std_success_durations = zeros(1,num_buckets)';
std_fail_durations = zeros(1,num_buckets)';
for i=1:num_buckets
    mean_success_durations(i) = mean(success_durations{i});
    mean_fail_durations(i) = mean(fail_durations{i});
    std_success_durations(i) = std(success_durations{i});
    std_fail_durations(i) = std(fail_durations{i});
end
mean_success_durations
std_success_durations

figure
plot(heights, mean_success_durations./24,'kx',heights,mean_fail_durations./24,'kx')
hold on
% Put some error bars on
for i=1:num_buckets
    % error bars
    plot([heights(i) heights(i)],[mean_success_durations(i)-std_success_durations(i) mean_success_durations(i)+std_success_durations(i)]./24,'k--')
    plot([heights(i) heights(i)],[mean_fail_durations(i)-std_fail_durations(i) mean_fail_durations(i)+std_fail_durations(i)]./24,'k-')
    % lines across top and bottom of bars!
    plot([heights(i)+.1 heights(i)-.1],[mean_fail_durations(i)-std_fail_durations(i) mean_fail_durations(i)-std_fail_durations(i)]./24,'k-')
    plot([heights(i)+.1 heights(i)-.1],[mean_fail_durations(i)+std_fail_durations(i) mean_fail_durations(i)+std_fail_durations(i)]./24,'k-')
    plot([heights(i)+.1 heights(i)-.1],[mean_success_durations(i)-std_success_durations(i) mean_success_durations(i)-std_success_durations(i)]./24,'k-')
    plot([heights(i)+.1 heights(i)-.1],[mean_success_durations(i)+std_success_durations(i) mean_success_durations(i)+std_success_durations(i)]./24,'k-')
end
set(gca,'FontSize',14)
xlim([-7 1])
ylim([0 100])
xlabel('log_{10} Initial mutation height (cell diameters)','FontSize',16)
ylabel('Mean simulation duration (days)','FontSize',16)

figure
bar(heights,fail_space_buckets,'k')
hold on
bar(heights,success_space_buckets,'r')
set(gca,'FontSize',14)
xlim([-7 1])
xlabel('log_{10} Initial mutation height (cell diameters)','FontSize',16)
ylabel('Number of experiments','FontSize',16)
legend('Failures','Dominations','Location','NorthWest')

% figure
% bar(heights,100*probability_space_buckets,'k')
% set(gca,'FontSize',14)
% xlabel('log10(Height of ancestor cell (cell diameters))','FontSize',16)
% ylabel('Probability of becoming dominant clone (%)','FontSize',16)
% 
% num_buckets = 30;
% success_times_in_days = data(success_indices,3)./24;
% x = linspace(0,max(success_times_in_days),num_buckets);
% n_elements = histc(success_times_in_days,x);
% c_elements = 100*cumsum(n_elements)./length(success_indices);
% 
% figure
% subplot(1,2,1)
% bar(x,n_elements,'k')
% set(gca,'FontSize',14)
% ylabel('Number of simulations','FontSize',16)
% xlabel('Time to monoclonal conversion (days)','FontSize',16)
% xlim([0 max(success_times_in_days)])
% subplot(1,2,2)
% plot(x,c_elements,'kx-')
% set(gca,'FontSize',14)
% xlim([0 max(success_times_in_days)])
% ylabel('Percentage of monoclonal crypts','FontSize',16)
% xlabel('Time (days)','FontSize',16)

