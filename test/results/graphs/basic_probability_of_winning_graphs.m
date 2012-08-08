function basic_probability_of_winning_graphs()

[data settings] = loadData();

% Load just the Wild-Type data (labelled no change to viscosity or
% proliferation) from all heights

% Locate the wild type data in the loaded data
wt_idx = -1;
for i=1:length(settings.viscosity_names)
    if strmatch(settings.viscosity_names{i},'1.0')
        wt_idx = i;
        break
    end
end

% Get data from all three heights of initial mutation.
data = [data{4,1,wt_idx}; data{4,2,wt_idx}; data{4,3,wt_idx}];

success_indices = find(data(:,4)==1);
failure_indices = find(data(:,4)==0);

length(data)

% Sanity check
assert(length(success_indices) + length(failure_indices)==length(data));

%num_successes = length(success_indices);

figure
subplot(2,1,1)
scatter(data(success_indices,2),data(success_indices,3)./24,'k.')
xlim([0 2.15])
ylim([0 250])
xlabel('Initial mutation height (cell diameters)','FontSize',16)
ylabel('Simulation duration (days)','FontSize',16)
set(gca,'FontSize',14)
subplot(2,1,2)
scatter(data(failure_indices,2),data(failure_indices,3)./24,'k.')
xlim([0 2.15])
ylim([0 250])
xlabel('Initial mutation height (cell diameters)','FontSize',16)
ylabel('Simulation duration (days)','FontSize',16)
set(gca,'FontSize',14)

num_buckets = 15;
heights = max(data(:,2))/num_buckets;
success_space_buckets = zeros(1,num_buckets);
fail_space_buckets = zeros(1,num_buckets);

for i = 1:length(data)
    for b = 1:num_buckets
        if ((b-1)*heights < data(i,2) && data(i,2) <= b*heights)
            if (data(i,4)==1)
                success_space_buckets(b) = success_space_buckets(b) + 1;
            else
                fail_space_buckets(b) = fail_space_buckets(b) + 1;
            end
        end
    end
end
probability_space_buckets = success_space_buckets ./...
    (success_space_buckets + fail_space_buckets);

figure
bar(([1:num_buckets]*heights),100*probability_space_buckets,'k')
set(gca,'FontSize',14)
xlabel('Height of ancestor cell (cell diameters)','FontSize',16)
ylabel('Probability of becoming dominant clone (%)','FontSize',16)

num_buckets = 30;
success_times_in_days = data(success_indices,3)./24;
failure_times_in_days = data(failure_indices,3)./24;
x = linspace(0,max([success_times_in_days; failure_times_in_days]),num_buckets);
n_elements = histc(success_times_in_days,x);
nf_elements = histc(failure_times_in_days,x);
c_elements = 100*cumsum(n_elements)./length(success_indices);

figure
subplot(1,2,1)
bar(x,n_elements,'k')
set(gca,'FontSize',14)
ylabel('Number of simulations','FontSize',16)
xlabel('Time to monoclonal conversion (days)','FontSize',16)
xlim([0 max(success_times_in_days)])
subplot(1,2,2)
plot(x,c_elements,'kx-')
set(gca,'FontSize',14)
xlim([0 max(success_times_in_days)])
ylabel('Percentage of monoclonal crypts','FontSize',16)
xlabel('Time (days)','FontSize',16)

greyness = .70;

% figure
% subplot(1,2,1)
% bar(x,nf_elements,0.7,'k')
% hold on
% bar(x-1,n_elements,0.7,'FaceColor',[greyness,greyness,greyness])
% set(gca,'FontSize',14)
% ylabel('Number of simulations','FontSize',16)
% xlabel('Simulation duration (days)','FontSize',16)
% xlim([-x(2) max(failure_times_in_days)])
% subplot(1,2,2)
% bar(x,nf_elements,0.7, 'k')
% hold on
% bar(x-1,n_elements,0.7,'FaceColor',[greyness,greyness,greyness])
% set(gca,'FontSize',14)
% ylabel('Number of simulations','FontSize',16)
% xlabel('Simulation duration (days)','FontSize',16)
% xlim([-x(2) max(failure_times_in_days)])
% ylim([0 120])
% legend('Clone lost','Clone dominates','Location','NorthEast')

% I think we really want a 3D bar chart here...

figure
bar3(x,[log10(n_elements) log10(nf_elements)],'Detached')
ylim([-10 160])
set(gca,'FontSize',14)
legend('Clone dominates','Clone lost','EastOutside')

ylabel('Simulation duration (days)','FontSize',16)
zlabel('Number of simulations','FontSize',16)


% 
% figure
% bar(x,log10(nf_elements),0.7,'k')
% %set(gca,'YScale','log')
% hold on
% bar(x,-log10(n_elements),0.7,'FaceColor',[greyness,greyness,greyness])
% set(gca,'FontSize',14)
% ylabel('Number of simulations','FontSize',16)
% xlabel('Simulation duration (days)','FontSize',16)
% xlim([-x(2) max(failure_times_in_days)])
% ylim([-4.1 4.1])
% legend('Clone lost','Clone dominates','Location','NorthEast')

% figure
% bar(x,nf_elements,0.7,'k')
% set(gca,'YScale','log')
% hold on
% bar(x-1,-n_elements,0.7,'FaceColor',[greyness,greyness,greyness])
% set(gca,'FontSize',14)
% ylabel('Number of simulations','FontSize',16)
% xlabel('Simulation duration (days)','FontSize',16)
% xlim([-x(2) max(failure_times_in_days)])
% ylim([-120 inf])
% legend('Clone lost','Clone dominates','Location','NorthEast')
