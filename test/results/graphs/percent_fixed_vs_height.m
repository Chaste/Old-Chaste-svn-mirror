function [ output_args ] = percent_fixed_vs_height( input_args )
%PERCENT_FIXED_VS_HEIGHT Summary of this function goes here
%   Detailed explanation goes here
[data, settings] = loadData();
num_heights = length(settings.heights);

num_buckets = 12;

% These are hard-coded into the simulations
crypt_length = 14.500;
box_proportion = 1.0/20.0;

for v=1:length(settings.viscosities)
    figure
    for i=1:length(settings.mutations)
        dataset=[];
        for j=1:length(settings.heights)
            % For all heights and mutations (viscosity 2 = 1.0, control)
            dataset = [dataset;data{i,j,v}];
        end
        data{i,v} = dataset;
        if ~isempty(dataset)
            takeover_index = find(dataset(:,4)==1);
            failure_index = find(dataset(:,4)==0);

            ordering = [2 4 3 1];
            subplot(2,2,ordering(i))
            scatter(dataset(failure_index,2), dataset(failure_index,3)./24,'bo')
            hold on
            scatter(dataset(takeover_index,2), dataset(takeover_index,3)./24,'rx')
            set(gca,'Fontsize',14)
            xlabel('Initial mutation height (cell diameters)','Fontsize',16)
            ylabel('Simulation duration (days)','Fontsize',16)
            title([settings.mutation_names{i} ', \alpha = ' settings.viscosity_names{v}],'Fontsize',16)
            if i==1
                legend('Mutation swept out','Mutation dominates','Location','NorthEast')
            end
        end
    end
end

% 
% % Work out proportions
% top_height = box_proportion*crypt_length*num_heights;
% bucket_height = top_height/num_buckets;
% this_bucket_height_min = linspace(0,top_height-bucket_height,num_buckets)';
% this_bucket_height_max = this_bucket_height_min + bucket_height;
% 
% for v = 1:length(settings.viscosities)
%     Success=[];
%     Total=[];
%     for m = 1:length(settings.mutations)
%         dataset = data{m,v};
%         successes=zeros(1,num_buckets);
%         failures=zeros(1,num_buckets);
%         for b=1:num_buckets
%             for i=1:length(dataset)
%                 if b == num_buckets % Took me ages to figure this one out!
%                     upper_lim = (dataset(i,2)<=this_bucket_height_max(b)+1e-4);
%                 else
%                     upper_lim = (dataset(i,2)<this_bucket_height_max(b));
%                 end
%                 if (dataset(i,2)>=this_bucket_height_min(b) && upper_lim )
%                     if dataset(i,4)==1
%                         successes(b) = successes(b)+1;
%                     else
%                         failures(b) = failures(b) +1;
%                     end
%                 end
%             end
%         end
%         assert(sum(successes)+sum(failures)==length(dataset));
%         
%         Success(:,m) = successes;
%         Total(:,m) = successes+failures;
%     end
%     
%     AllMutationsSuccess = sum(Success,2);
%     AllMutationsTotal = sum(Total,2);
%     
%     Proportions_all_mutations = 100.0*AllMutationsSuccess./AllMutationsTotal;
%     
%     for m=1:length(settings.mutations)
%         Proportion_of_Success_from_each_mutation(:,m) = Proportions_all_mutations.*Success(:,m)./sum(Success,2);
%     end
%     
%     figure
%     subplot(1,2,1)
%     bar(this_bucket_height_max,Proportions_all_mutations)
%     set(gca,'FontSize',14)
%     xlabel('Height below which cell was labelled','FontSize',16)
%     ylabel('Percentage which will dominate the crypt','FontSize',16)
%     title(['Origin of dominating cell, by mutation height, viscosity = ' settings.viscosity_names{v}],'FontSize',16)
% 
%     Proportions_for_all_viscosities(:,v) = Proportions_all_mutations;
%     
%     subplot(1,2,2)
%     bar(this_bucket_height_max,Proportion_of_Success_from_each_mutation,'Stack')
%     set(gca,'FontSize',14)
%     xlabel('Height below which cell was labelled','FontSize',16)
%     ylabel('Percentage which will dominate the crypt','FontSize',16)
%     title('Which mutations contrbute to this','FontSize',16)
%     legend(settings.mutation_names)
%     
%     if v==3
%         figure(123)
%         subplot(1,2,1)
%         bar(this_bucket_height_max,100.0*Success./Total,'Group')
%         set(gca,'FontSize',14,'XTick',round(100*this_bucket_height_max)/100.0)
%         xlim([0 1.5])
%         xlabel('Height below which cell was labelled','FontSize',16)
%         ylabel('Percentage which dominate the crypt','FontSize',16)
%         title('Effect of varying proliferation of mutation','FontSize',16)
%         legend(settings.mutation_names)
%     end
%     
% end
% 
% figure(123)
% subplot(1,2,2)
% bar(this_bucket_height_max,Proportions_for_all_viscosities,'Group')
% set(gca,'FontSize',14,'XTick',round(100*this_bucket_height_max)/100.0)
% xlim([0 1.5])
% xlabel('Height below which cell was labelled','FontSize',16)
% ylabel('Percentage which dominate','FontSize',16)
% title('Effect of varying viscosity of mutation','FontSize',16)
% legend(settings.viscosity_names)
% 
% %% Export this figure to the paper images directory.
% papersize = get(gcf, 'PaperSize');
% width = 30.0;         % Initialize a variable for width.
% height = 12.0;      % Initialize a variable for height.
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(gcf, 'PaperPosition', myfiguresize);
% 
% print -depsc2 ../../../paper/images/proliferation_vs_viscosity



