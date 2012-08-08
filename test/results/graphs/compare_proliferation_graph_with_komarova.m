function [ Z ] = compare_proliferation_graph_with_komarova( input_args )


% Get the 'control' simulation from old results
data = loadData();
control_data = data{4,1,7};
control_prob = sum(control_data(:,4))./length(control_data);

%   Detailed explanation goes here
[data, settings] = load_proliferation_only_data();

% Put the height into sub-boxes...
num_heights = length(settings.proliferation_heights);

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

% Plot of probabilities for different heights
success_prob = (success_bucket./num_bucket);
% Hack in the control case
success_prob = [success_prob(1); control_prob; success_prob(2:end)];

%confidence_interval = 0.975;
%errors = norminv(confidence_interval,0,1) .* sqrt((success_prob.*(1-success_prob))./num_bucket);

N = 16;
r = linspace(0,3.5,10000);
rho = 100* ( 2*(r.^(N-1)).*(1-r)./(1 + r + r.^(N-1) - 3*r.^N) );

% Plot the expression for the probability of mutant domination on its own
figure
hold on
plot(r,rho,'k-','LineWidth',1)
xlabel('Mutant fitness{\it r}','FontSize',16)
ylabel('Probability (%)','FontSize',16)
set(gca,'FontSize',16)
axis([0 3.5 0 60])
box on

success_percentage = 100*success_prob;
r_values = zeros(size(success_percentage));
for i=1:length(success_percentage)
    % Find the value of the fitness r at which the domination probabilitiy
    % takes this value
    p = success_percentage(i);
    idx = find(rho>p);    
    r_values(i) = r(idx(1));
    %r_values(i) = fsolve(@(r) p-100*(2*(r.^13).*(1-r)./(1+r+ r.^13 - 3*r.^14) ),[0.1],optimset('Display','off'))
end
alpha_values = this_bucket_height;

% Plot the alpha values against the 'corresponding' r values to check if
% there is a simple relationship between these parameters
figure
hold on
heights = 100*[this_bucket_height(1); 0.35; this_bucket_height(2:end)];
plot(heights, r_values,'ko','LineWidth',1)
ylabel('Mutant fitness{\it r}','FontSize',16)
xlabel('Proliferation ceiling{\it y}_{thr} (% of crypt)','FontSize',16)
set(gca,'FontSize',16)
% axis([0 11 0 70])
box on
p = polyfit(heights,r_values,1)
f = polyval(p,[25; heights]);
plot([25; heights], f,'k-','LineWidth',1)
ylim([0.6 1.3])
xlim([25 100])
