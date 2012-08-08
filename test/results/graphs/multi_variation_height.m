function [ output_args ] = multi_variation_height( input_args )
%PERCENT_FIXED_VS_HEIGHT Summary of this function goes here
%   Detailed explanation goes here
[data, settings] = loadData();

alpha_1_idx = -1;
alpha_2_idx = -1;
% Locate the alpha = 1 and alpha = 2 data in the loaded data
for i=1:length(settings.viscosity_names)
    if strmatch(settings.viscosity_names{i},'1.0')
        alpha_1_idx = i;
    end
        if strmatch(settings.viscosity_names{i},'2.0')
        alpha_2_idx = i;
    end
end

% These are hard-coded into the simulations
crypt_length = 14.500;
box_proportion = 1.0/20.0;
ordering = [2 4 3 1];
settings.mutations
for v=1:length(settings.viscosities)
   % figure
    for i=1:length(settings.mutations)
        dataset=[];
        for j=1:length(settings.heights)
            % For all heights and mutations (viscosity 7 = 1.0, control)
            dataset = [dataset;data{i,j,v}];
        end
        data{i,v} = dataset;
        
        % Do all the scatter plots
%         if ~isempty(dataset)
%             takeover_index = find(dataset(:,4)==1);
%             failure_index = find(dataset(:,4)==0);
% 
%             subplot(2,2,ordering(i))
%             scatter(dataset(failure_index,2), dataset(failure_index,3)./24,'bo')
%             hold on
%             scatter(dataset(takeover_index,2), dataset(takeover_index,3)./24,'rx')
%             set(gca,'Fontsize',14)
%             xlabel('Initial mutation height (cell diameters)','Fontsize',16)
%             ylabel('Simulation duration (days)','Fontsize',16)
%             title([settings.mutation_names{i} ', \alpha = ' settings.viscosity_names{v}],'Fontsize',16)
%             if i==1
%                 legend('Mutation swept out','Mutation dominates','Location','NorthEast')
%             end
%         end
        
    end
end

% These correspond to alpha = 1 and 2 as examples.
max_v = alpha_2_idx;
min_v = alpha_1_idx;

figure
for v=[min_v max_v]
    for m=1:length(settings.mutations)
        this_data = data{m,v};
        length(this_data)
        num_buckets = 15;
        num_duration_buckets = 30;
        heights = max(this_data(:,2))/num_buckets;
        durations = max(this_data(:,3))/num_duration_buckets;
        success_space_buckets = zeros(1,num_buckets);
        fail_space_buckets = zeros(1,num_buckets);
        success_duration_buckets = zeros(1,num_duration_buckets);
        fail_duration_buckets = zeros(1,num_duration_buckets);
        
        
        for i = 1:length(this_data)
            for b = 1:num_buckets
                if ((b-1)*heights < this_data(i,2) && this_data(i,2) <= b*heights)
                    if (this_data(i,4)==1)
                        success_space_buckets(b) = success_space_buckets(b) + 1;
                    else
                        fail_space_buckets(b) = fail_space_buckets(b) + 1;
                    end
                end
            end
            
            for b = 1:num_duration_buckets
                if ((b-1)*durations < this_data(i,3) && this_data(i,3) <= b*durations)
                    if (this_data(i,4)==1)
                        success_duration_buckets(b) = success_duration_buckets(b) + 1;
                    else
                        fail_duration_buckets(b) = fail_duration_buckets(b) + 1;
                    end
                end
            end
        end       
        
        %% Probability vs height
        probability_space_buckets = success_space_buckets ./ (success_space_buckets + fail_space_buckets);
        if v==max_v
            xplot1(:,m,2) = [1:num_buckets]*heights;
            yplot1(:,m,2) = 100*probability_space_buckets;
            xplot2(:,m,2) = [1:num_duration_buckets]*durations./24;
            yplot2(:,m,2) = success_duration_buckets./length(find(this_data(:,4)==1));
            
        else
            xplot1(:,m,1) = [1:num_buckets]*heights;
            yplot1(:,m,1) = 100*probability_space_buckets;
            xplot2(:,m,1) = [1:num_duration_buckets]*durations./24;
            yplot2(:,m,1) = success_duration_buckets./length(find(this_data(:,4)==1));        
        end
        
%         figure(1)
%         subplot(2,4,ordering(m)+offset)
%         bar([1:num_buckets]*heights, 100*probability_space_buckets,'k')
%         set(gca,'FontSize',12)
%         ylim([0 30])
%         xlim([0 2])
%         if v==max_v
%         xlabel('Initial mutation height (cell diameters)','FontSize',14)
%         end
%         if m==4
%         ylabel('Probability (%)','FontSize',14)
%         end
%         title([settings.mutation_names{m} ', \alpha = ' settings.viscosity_names{v}],'Fontsize',14)
%         %% Number running to what duration
%         
%         figure(2)
%         subplot(2,4,ordering(m)+offset)
%         bar([1:num_duration_buckets]*durations./24,success_duration_buckets./length(find(this_data(:,4)==1)),'k')
%         xlim([0 8000./24])
%         ylim([0 0.25])
%         if v==max_v
%             xlabel('Duration of dominant clones (days)','FontSize',14)
%         end
%         if m==4
%             ylabel('Proportion of simulations','FontSize',14)
%         end
%         title([settings.mutation_names{m} ', \alpha = ' settings.viscosity_names{v}],'Fontsize',14)
    end
end

figure(1)
createfigure_multi_variation(xplot1, yplot1,'Mutation height (cells)', 'Probability (%)',2,31);

figure(2)
createfigure_multi_variation(xplot2, 100.*yplot2,{'Time for clone to become','dominant (days)'}, '% of simulations',250,25);

