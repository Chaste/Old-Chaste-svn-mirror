
function [ Z ] = compare_adhesion_graph_with_komarova( input_args )

% This function reproduces much of the function probability_of_domination_adhesion_only
[data, settings] = loadData();

for i=1:length(settings.viscosity_names)
    data_reduced{i} = data{4,1,i};
end

num_viscos = length(settings.viscosity_names);
num_buckets = 1;

this_bucket_height = [];
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
            i = j + (viscosity_idx-1)*num_buckets;
            this_bucket_dataset = dataset;
            num_bucket(i) = num_bucket(i) + length(dataset);
            success_bucket(i) = success_bucket(i) + length(find(this_bucket_dataset(:,4)==1));
            failure_bucket(i) = failure_bucket(i) + length(find(this_bucket_dataset(:,4)==0));
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

success_prob = (success_bucket./num_bucket);
failure_prob = (failure_bucket./num_bucket);

confidence_interval = 0.975;
errors = norminv(confidence_interval,0,1) .* sqrt((success_prob.*(1-success_prob))./num_bucket);

N = 16;
r = linspace(0,3.5,100000);
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
plot(alpha_values, r_values,'ko','LineWidth',1)
ylabel('Mutant fitness{\it r}','FontSize',16)
xlabel('Mutant adhesion parameter \alpha','FontSize',16)
set(gca,'FontSize',16)
% axis([0 11 0 70])
box on
p = polyfit(alpha_values,r_values,1)
f = polyval(p,[0; alpha_values]);
plot([0; alpha_values], f,'k-','LineWidth',1)
