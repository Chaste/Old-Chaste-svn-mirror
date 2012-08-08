function basic_wild_type_graphs()

data = load_WT_Data();

% Does the original cell ever remain when crypt is monoclonal?
well = find(data(:,3)>data(:,2));
% No.
assert(isempty(well))

num_buckets = 40;
success_times_in_days = data(:,2)./24;
x = linspace(0,max(success_times_in_days)+10,num_buckets);
n_elements = histc(success_times_in_days,x);
c_elements = 100*cumsum(n_elements)./length(success_times_in_days);

% What is longest time taken?
max_time_taken = (max(data(:,2))./24)./7

figure
subplot(1,2,1)
bar(x,n_elements,'k')
set(gca,'FontSize',14)
ylabel('Number of simulations','FontSize',16)
xlabel('Time to monoclonal conversion (days)','FontSize',16)
xlim([0 max(x)])
subplot(1,2,2)
plot(x,c_elements,'kx-')
set(gca,'FontSize',14)
xlim([0 max(x)])
ylabel('Monoclonal crypts (%)','FontSize',16)
xlabel('Time (days)','FontSize',16)


num_buckets = 40;
success_times_in_days = data(:,3)./24;
x = linspace(0,max(success_times_in_days)+2,num_buckets);
n_elements = histc(success_times_in_days,x);
c_elements = 100*(1.0-cumsum(n_elements)./length(success_times_in_days));

figure
subplot(1,2,1)
bar(x,n_elements,'k')
set(gca,'FontSize',14)
ylabel('Number of simulations','FontSize',16)
xlabel('Time to ancestor cell loss (days)','FontSize',16)
xlim([0 max(x)])
subplot(1,2,2)
plot(x,c_elements,'kx-')
set(gca,'FontSize',14)
xlim([0 max(x)])
ylabel('Crypts retaining ancestor (%)','FontSize',16)
xlabel('Time (days)','FontSize',16)