function [ output_args ] = time_to_domination( input_args )

% Time-to-monoclonality plots
[data, settings] = loadData();

type_dist=[];
type_dist_2=[];
% Pull out the monoclonal conversion time
% (dominate only count)
for state=1:length(settings.mutations)
    for height = 1:length(settings.heights)
        for visc = 1:length(settings.viscosities)
            try
                dataset = data{state,height,visc};
                
                if visc==3
                    type_dist = [type_dist; [dataset state.*ones(1,length(dataset))']];
                end
            
                type_dist_2 = [type_dist_2; [dataset visc.*ones(1,length(dataset))' state.*ones(1,length(dataset))' ]];
            catch
                disp('Data file doesn''t exist') 
            end

        end
    end
end
dominate = find(type_dist(:,4)==1);
type_dist = type_dist(dominate,:);
% length(type_dist)

dominate = find(type_dist_2(:,4)==1);
type_dist_2 = type_dist_2(dominate,:);

top_time = ceil(max([type_dist(:,3); type_dist_2(:,3)]));

buckets = zeros(length(settings.mutations),top_time);
cumulative_buckets = buckets;

buckets_1 = zeros(length(settings.mutations),length(settings.viscosities),top_time);
cumulative_buckets_1 = buckets_1;

% size(buckets)
% size(buckets_1)

for i=1:length(buckets)
    
    for j=1:length(type_dist)
        if type_dist(j,3) > i-1 && type_dist(j,3) <= i
            if type_dist(j,4)==1 % If this is a domination
                buckets(type_dist(j,5),i) = buckets(type_dist(j,5),i) +1 ;
            end
        end
    end
    
    for j=1:length(type_dist_2)
        if type_dist_2(j,3) > i-1 && type_dist_2(j,3) <= i
            if type_dist_2(j,4)==1 % If this is a domination
                buckets_1(type_dist_2(j,6),type_dist_2(j,5),i) = buckets_1(type_dist_2(j,6),type_dist_2(j,5),i) +1 ;
            end
        end
    end
    
    if i>1
        cumulative_buckets(:,i) = buckets(:,i) + cumulative_buckets(:,i-1);
        for j=1:length(settings.mutations)
            cumulative_buckets_1(j,:,i) = buckets_1(j,:,i) + cumulative_buckets_1(j,:,i-1);
        end
    else
        cumulative_buckets(:,i) = buckets(:,i);
        for j=1:length(settings.mutations)
            cumulative_buckets_1(j,:,i) = buckets_1(j,:,i);
        end
    end
end

figure
% Need to think of a nice way to present this in 3D maybe...
% Just the wild-type
colours = {'r-','b-','k:','g-','m-','r--','k-','b--','g--','k--','m--','g:','m:','r:'};
for v=1:length(settings.viscosities)
    plot_this = 100.0*cumulative_buckets_1(4,v,:)./cumulative_buckets_1(4,v,end);
    plot_this = squeeze(plot_this);
    plot3(1:top_time, settings.viscosities(v).*ones(top_time,1),plot_this,colours{v})
    hold on
end
zlim([0 100])
xlim([0 top_time])
xlabel('Time to domination (hours)','FontSize',16)
zlabel('Cumulative percentage of simulations','FontSize',16)
ylabel('Viscosity \alpha','FontSize',16)
set(gca,'FontSize',14)
title(['Mutation State = ' settings.mutation_names{4}],'FontSize',16)


%% FIRST FIGURE
figure
for i=1:length(settings.mutations)
    subplot(2,2,i)
    colours = {'r-','b-','k:','g-','m-','r--','k-','b--','g--','k--','m--','g:','m:','r:'};
    for v=1:length(settings.viscosities)
        plot_this = 100*cumulative_buckets_1(i,v,:)./cumulative_buckets_1(i,v,end);
        plot_this = squeeze(plot_this);
        plot(1:top_time,plot_this,colours{v})
        hold all
    end
    xlim([0 top_time])
    ylim([0 100])
    xlabel('Time to domination (hours)','FontSize',16)
    ylabel('Cumulative percentage of simulations','FontSize',16)
    set(gca,'FontSize',14)
    legend(settings.viscosity_names,'Location','SouthEast')
    title(['Mutation State = ' settings.mutation_names{i}],'FontSize',16)
end
%% Export this figure to the paper images directory.
papersize = get(gcf, 'PaperSize');
width = 30.0;         % Initialize a variable for width.
height = 22.0;      % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

print -depsc2 ../../../paper/images/time_to_domination_all_WT


%% SECOND FIGURE
figure
colours = {'r--','g-','b-.','k:'};
for m=1:4
    plot(1:top_time,100*cumulative_buckets(m,:)./cumulative_buckets(m,end),colours{m})
    hold all
end
xlim([0 top_time])
ylim([0 100])
xlabel('Time to domination of crypt (hours)','FontSize',16)
ylabel('Cumulative percentage of simulations','FontSize',16)
set(gca,'FontSize',14)
legend(settings.mutation_names,'Location','SouthEast')

%% Export this figure to the paper images directory.
papersize = get(gcf, 'PaperSize');
width = 15.0;         % Initialize a variable for width.
height = 10.0;      % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

print -depsc2 ../../../paper/images/time_to_domination_same_viscos

